# Knuth-Morris-Pratt (KMP) Algorithm
def kmp_search(text, pattern):
    def compute_lps(pattern):
        lps = [0] * len(pattern)
        length = 0
        i = 1
        while i < len(pattern):
            if pattern[i] == pattern[length]:
                length += 1
                lps[i] = length
                i += 1
            else:
                if length != 0:
                    length = lps[length - 1]
                else:
                    lps[i] = 0
                    i += 1
        return lps

    lps = compute_lps(pattern)
    i = 0  # index for text
    j = 0  # index for pattern
    positions = []
    
    while i < len(text):
        if pattern[j] == text[i]:
            i += 1
            j += 1

        if j == len(pattern):
            positions.append(i - j)
            j = lps[j - 1]
        elif i < len(text) and pattern[j] != text[i]:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1
    return positions

# Rabin-Karp Algorithm
def rabin_karp_search(text, pattern, d=256, q=101):
    n = len(text)
    m = len(pattern)
    h = 1
    p = 0  # hash value for pattern
    t = 0  # hash value for text
    positions = []

    for i in range(m-1):
        h = (h * d) % q

    for i in range(m):
        p = (d * p + ord(pattern[i])) % q
        t = (d * t + ord(text[i])) % q

    for i in range(n - m + 1):
        if p == t:
            if text[i:i + m] == pattern:
                positions.append(i)
        if i < n - m:
            t = (d * (t - ord(text[i]) * h) + ord(text[i + m])) % q
            if t < 0:
                t = t + q
    return positions

# Boyer-Moore Algorithm
def boyer_moore_search(text, pattern):
    def bad_character_heuristic(pattern):
        bad_char = [-1] * 256
        for i in range(len(pattern)):
            bad_char[ord(pattern[i])] = i
        return bad_char

    m = len(pattern)
    n = len(text)
    bad_char = bad_character_heuristic(pattern)
    s = 0
    positions = []

    while s <= n - m:
        j = m - 1

        while j >= 0 and pattern[j] == text[s + j]:
            j -= 1

        if j < 0:
            positions.append(s)
            s += (m - bad_char[ord(text[s + m])] if s + m < n else 1)
        else:
            s += max(1, j - bad_char[ord(text[s + j])])
    
    return positions

# Example usage:
text = "ABABDABACDABABCABAB"
pattern = "ABABCABAB"

print("KMP Search:", kmp_search(text, pattern))
print("Rabin-Karp Search:", rabin_karp_search(text, pattern))
print("Boyer-Moore Search:", boyer_moore_search(text, pattern))
