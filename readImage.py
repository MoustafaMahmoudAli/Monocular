from PIL import Image

im = Image.open("F:\\ENG\\MONO ROBO\\test\\samples\\horse.png")
im = im.convert("RGB")
w, h = im.size
f = open("frame1.txt", "w")
f.write(str(w))
f.write(" ")
f.write(str(h))
f.write(" ")
for j in range(0,h):
    for i in range(0,w):
        p = im.getpixel((i,j))
        f.write(str(p[0]))
        f.write(" ")
        f.write(str(p[1]))
        f.write(" ")
        f.write(str(p[2]))
        f.write(" ")
f.close()

