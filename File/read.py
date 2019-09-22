import sys

def main():
    argv = sys.argv
    argc = len(argv)

    all_data = []
    for i in range(1,argc):
        with open(argv[i],mode = "r") as fp:
            lines = fp.readlines()
            for i in range(len(lines)):
                lines[i] = lines[i].strip()
            all_data.append(lines)

    with open ("all_data" , mode = "w") as fw:
        fw.write("    name  |  old\n")
        for i in range (len (all_data)):
            for j in range (len(all_data[i])):
                fw.write("{:>8}".format(all_data[i][j]))
            fw.write("\n")

if __name__=="__main__":
    main()
