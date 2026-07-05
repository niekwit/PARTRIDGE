from itertools import groupby

CIGAR_CONSUMES_QUERY = ("M", "I", "S", "=", "X", "L", "H")
CIGAR_CONSUMES_REFERENCE = ("M", "D", "N", "=", "X", "L")
CIGAR_CONSUMES_BOTH = ("M", "=", "X", "L")
CIGAR_QUERY_MATCHES = ("M", "=", "L")


class CigarItem:
    """
    Custom Definitions
    L = Maps to LTR (equal to M)
    """

    def __init__(self, length: int, operator: chr):
        self.operator = operator
        self.length = length
        self.query_offset = None
        self.reference_offset = None

    def __len__(self):
        return self.length

    def offsetIterFunc(self, offsets: (int, int)) -> (int, int):
        self.query_offset, self.reference_offset = offsets
        return self.query_offset + self.query_len, self.reference_offset + self.ref_len

    def __str__(self):
        return f"{self.length}{self.operator}"

    @property
    def ref_len(self):
        if self.operator in CIGAR_CONSUMES_REFERENCE:
            return self.length
        else:
            return 0

    @property
    def query_len(self):
        if self.operator in CIGAR_CONSUMES_QUERY:
            return self.length
        else:
            return 0

    def expand(self) -> str:
        return "".join([self.operator] * self.length)

    def __mul__(self, other: "CigarItem") -> int:
        """
        Calculate Query Overlap
        :param other:
        :return:
        """
        if (
            self.operator in CIGAR_QUERY_MATCHES
            and other.operator in CIGAR_QUERY_MATCHES
        ):
            if self.query_offset == other.query_offset:
                return min(self.query_len, other.query_len)
            elif self.query_offset > other.query_offset:
                if (
                    self.query_offset + self.query_len
                    >= other.query_offset + other.query_len
                ):
                    return max(
                        0, other.query_offset + other.query_len - self.query_offset
                    )
                else:
                    return self.query_len
            else:
                if (
                    other.query_offset + other.query_len
                    >= self.query_offset + self.query_len
                ):
                    return max(
                        0, self.query_offset + self.query_len - other.query_offset
                    )
                else:
                    return other.query_len
        return 0


class Cigar:
    def __init__(self, cigar_string: str, reverse=False):
        if cigar_string is None:
            self.cigar = []
            return
        cigar = groupby(cigar_string, lambda c: c.isdigit())
        self.cigar = []
        l = None
        for g, n in cigar:
            if l is None:
                l = int("".join(n))
                continue
            else:
                # if reverse:
                #    self.cigar.insert(0, CigarItem(l, ("".join(n))[0]))
                # else:
                self.cigar.append(CigarItem(l, ("".join(n))[0]))
                l = None
        self.setOffsets()

    def startend(self) -> (int, int):
        """
        Returns a tuple of integers indicating 0-based start and end-position of the part of the query matching the
        reference, relative to the query. These coordinates are always relative to the cigar string.
        :return:
        """
        q = 0
        left = None
        right = 0
        for c in self.cigar:
            if c.operator in CIGAR_CONSUMES_BOTH:
                if left is None:
                    left = q
                q += c.length
                right = q
            elif c.operator in CIGAR_CONSUMES_QUERY:
                q += c.length
        if left is None:
            return -1, 0
        else:
            return left, right

    def clipped(self) -> int:
        """
        Returns the number of soft-clipped bases
        :return: number of soft-clipped bases
        """
        clipped = 0
        for c in self.cigar:
            if c.operator == "S":
                clipped += c.length
        return clipped

    def purge(self):
        """
        Join adjacent identical operators together.
        :return:
        """
        last_c = None
        new_cigar = []
        for c in self.cigar:
            if last_c is None or last_c.operator != c.operator:
                if last_c is not None:
                    new_cigar.append(last_c)
                last_c = c
            else:
                last_c.length += c.length
        if last_c is not None:
            new_cigar.append(last_c)
        self.cigar = new_cigar

    def __str__(self):
        self.purge()
        return "".join([str(c) for c in self.cigar])

    def __eq__(self, other: "Cigar") -> bool:
        return str(self) == str(other)

    def setOffsets(self):
        offset = (0, 0)
        for c in self.cigar:
            offset = c.offsetIterFunc(offset)

    @property
    def ref_len(self) -> int:
        len = 0
        for c in self.cigar:
            len += c.ref_len
        return len

    @property
    def query_len(self) -> int:
        len = 0
        for c in self.cigar:
            len += c.query_len
        return len

    def expand(self) -> str:
        return "".join([c.expand() for c in self.cigar])

    @classmethod
    def collapse(cls, cigar: str):
        r = cls("")
        for c, n in groupby(cigar):
            n = list(n)
            r.cigar.append(CigarItem(len(n), c))
        r.setOffsets()
        return r

    def append(self, other: "Cigar"):
        self.cigar = self.cigar + other.cigar

    def __mul__(self, other: "Cigar") -> int:
        """
        Calculate query overlap
        :param other:
        :return:
        """
        r = 0
        for s in self.cigar:
            for o in other.cigar:
                r += s * o
        return r

    def is_disjoint(self, other: "Cigar", max_overlap=0) -> int:
        assert (
            self.query_len == other.query_len
        ), "Query lengths of self and other have to be equal, but they are not: {self.query_len} is not {other.query_len}"
        return self * other <= max_overlap

    def mark_ltr(self, offset: int):
        i = 0
        p = 0
        if offset >= 0:
            while i < len(self.cigar):
                if self.cigar[i].operator != "M":
                    p += self.cigar[i].ref_len
                    i += 1
                    continue
                if p > offset:
                    self.cigar[i].operator = "L"
                    i += 1
                    continue
                if self.cigar[i].ref_len + p <= offset:
                    i += 1
                    p += self.cigar[i].ref_len
                    continue
                not_matched = offset - p
                matched = self.cigar[i].ref_len - not_matched
                assert matched + not_matched == self.cigar[i].length
                if not_matched > 0:
                    self.cigar.insert(i, CigarItem(not_matched, self.cigar[i].operator))
                    self.cigar[i + 1].length = matched
                    self.cigar[i + 1].operator = "L"
                    p += matched + not_matched
                    i += 2
                else:
                    self.cigar[i].operator = "L"
                    p += self.cigar[i].ref_len
                    i += 1
        else:
            offset = self.ref_len + offset
            while i < len(self.cigar):
                if self.cigar[i].operator != "M":
                    p += self.cigar[i].ref_len
                    i += 1
                    continue
                print(p, self.cigar[i].ref_len, offset)
                if p + self.cigar[i].ref_len < offset:
                    self.cigar[i].operator = "L"
                    p += self.cigar[i].ref_len
                    i += 1
                    continue
                if p >= offset:
                    break
                not_matched = offset - p
                matched = self.cigar[i].ref_len - not_matched
                assert matched + not_matched == self.cigar[i].length
                if matched > 0:
                    self.cigar.insert(i, CigarItem(matched, "L"))
                    self.cigar[i + 1].length = not_matched
                    p += matched + not_matched
                    i += 2
                else:
                    self.cigar[i].operator = "L"
                    p += self.cigar[i].ref_len
                    i += 1
        self.setOffsets()

    def prefix(self, n):
        if n > 0:
            self.cigar.insert(0, CigarItem(n, "S"))
        return self

    def postfix(self, n):
        if n > 0:
            self.cigar.append(CigarItem(n, "S"))
        return self

    def add(self, operator: chr, length: int):
        assert length > 0, f"attempting to add cigar of length 0, operator={operator}"
        self.cigar.append(CigarItem(length, operator))

    def __add__(self, other: "Cigar"):
        """
        Merge two Cigars
        :param other: other cigar string
        :return:
        """
        assert (
            self.query_len == other.query_len
        ), f"Query length of cigar1 {self} does not correspond to query length of cigar 2 {other}"
        cigar1 = self.expand()
        cigar2 = other.expand()
        summary = ""
        i1 = 0
        i2 = 0
        q1 = 0
        q2 = 0
        if cigar2[i2] in CIGAR_CONSUMES_QUERY:
            q2 += 1
        if cigar1[i1] in CIGAR_CONSUMES_QUERY:
            q1 += 1
        while i1 < len(cigar1) and i2 < len(cigar2):
            if q1 > q2:
                summary += cigar2[i2]
                i2 += 1
                if i2 < len(cigar2) and cigar2[i2] in CIGAR_CONSUMES_QUERY:
                    q2 += 1
                continue
            if q2 > q1:
                summary += cigar1[i1]
                i1 += 1
                if i1 < len(cigar1) and cigar1[i1] in CIGAR_CONSUMES_QUERY:
                    q1 += 1
                continue
            if cigar1[i1] == cigar2[i2]:
                summary += cigar1[i1]
                i1 += 1
                i2 += 1
                if i2 < len(cigar2) and cigar2[i2] in CIGAR_CONSUMES_QUERY:
                    q2 += 1
                if i1 < len(cigar1) and cigar1[i1] in CIGAR_CONSUMES_QUERY:
                    q1 += 1
                continue
            if "L" == cigar1[i1] or "L" == cigar2[i2]:
                summary += "L"
                i1 += 1
                i2 += 1
                if i2 < len(cigar2) and cigar2[i2] in CIGAR_CONSUMES_QUERY:
                    q2 += 1
                if i1 < len(cigar1) and cigar1[i1] in CIGAR_CONSUMES_QUERY:
                    q1 += 1
                continue
            if "M" == cigar1[i1] or "M" == cigar2[i2]:
                summary += "M"
                i1 += 1
                i2 += 1
                if i2 < len(cigar2) and cigar2[i2] in CIGAR_CONSUMES_QUERY:
                    q2 += 1
                if i1 < len(cigar1) and cigar1[i1] in CIGAR_CONSUMES_QUERY:
                    q1 += 1
                continue
            summary += cigar1[i1]
            i1 += 1
            i2 += 1
            if i2 < len(cigar2) and cigar2[i2] in CIGAR_CONSUMES_QUERY:
                q2 += 1
            if i1 < len(cigar1) and cigar1[i1] in CIGAR_CONSUMES_QUERY:
                q1 += 1
        return Cigar.collapse(summary)

    def iterate(self, ref_pos: int = 0) -> (int, int, chr):
        """
        Iterate through a Cigar object.
        :param ref_pos: Reference position of the left most Cigar (note: By convention, cigars are ALWAYS relative to the + strand in these scripts, thus for a - alignment, set this to the last ref position, otherwise the first ref position)
        :param reverse_ref: set True if the read is opposite to the + strand (thus the reference is on the minus strand)
        :yield: A 3-tuple with query position (0-based), reference-position, and cigar char.
        """
        query_pos = 0
        for c in self.cigar:
            for i in range(c.length):
                yield ref_pos, query_pos, c.operator
                if c.operator in CIGAR_CONSUMES_QUERY:
                    query_pos += 1
                if c.operator in CIGAR_CONSUMES_REFERENCE:
                    ref_pos += 1

    def quickIterate(self, ref_pos: int = 0) -> (int, int, chr, int):
        """
        Same as iterate, but returns a forth parameter: actual length
        :param ref_pos:
        :return:
        """
        query_pos = 0
        for c in self.cigar:
            yield ref_pos, query_pos, c.operator, c.length
            if c.operator in CIGAR_CONSUMES_QUERY:
                query_pos += c.length
            if c.operator in CIGAR_CONSUMES_REFERENCE:
                ref_pos += c.length


if __name__ == "__main__":
    cigar1 = Cigar("16S34M5S")
    cigar2 = Cigar("16M39S")
    cigar3 = Cigar("18M37S")
    print(cigar1)
    assert cigar1 != cigar2
    assert cigar1 == Cigar.collapse(cigar1.expand())
    assert cigar2 == Cigar.collapse(cigar2.expand())
    assert cigar2 * cigar1 == cigar2 * cigar1
    assert cigar1.is_disjoint(cigar2)
    assert not cigar1.is_disjoint(cigar1)
    assert not cigar2.is_disjoint(cigar2)
    assert cigar1 * cigar3 == 2
    assert cigar1.is_disjoint(cigar2, 2)
    cigar1 = Cigar("16S34M2I4M5S")
    print(cigar1)
    cigar1.mark_ltr(2)
    print(cigar1)
    cigar1 = Cigar("16S34M2I4M5S")
    cigar1.mark_ltr(-2)
    print(cigar1)
