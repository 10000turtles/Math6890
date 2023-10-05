f = open("broke.txt", "r")

s = f.read()

spl = str.split(s, " ")

new_str = ""

key_words = ["static", "inline", "using", "double", "define",
             "ifndef", "include", "typedef", "const", "return"]

count = 2

for i in spl:
    try:
        if (int(i) == count):
            new_str = new_str + "\n"
            count = count + 1
    except:
        out = i
        for word in key_words:
            if word in out:
                index = i.find(word)+len(word)
                out = out[:index] + " " + out[index:]
        new_str = new_str + out

wr = open("poisson.C", "w")

wr.write(new_str)
