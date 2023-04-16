def get_nextval(patt):
    i = 0
    j = -1
    nextval = [-1]

    while i < len(patt)-1:
        if j == -1 or patt[i] == patt[j]:
            i += 1
            j += 1
            if patt[i] != patt[j]:
                nextval.append(j)
            else:
                nextval.append(nextval[j])
        else:
            j = nextval[j]
    return nextval

def kmp(text,patt):
    i = 0
    j = 0
    nextval = get_nextval(patt)

    while i < len(text) and j < len(patt):
        if j == -1 or text[i] == patt[j]:
            i += 1
            j += 1
        else:
            j = nextval[j]
    if j == len(patt):
        return i - j
    else:
        return -1



# 测试
# if __name__ == "__main__":
#     text = "aaacaaaabe"
#     pattern = "aaaab"

#     pos = kmp(text,pattern)
#     print(pos)