import struct

def main():
    with open("smoke.vol", "rb") as infile:
        data = infile.read()

        print(data[2])

        print(struct.unpack("<f", data[24:28]))
        print(struct.unpack("<f", data[28:32]))
        print(struct.unpack("<f", data[32:36]))
        print(struct.unpack("<f", data[36:40]))
        print(struct.unpack("<f", data[40:44]))
        print(struct.unpack("<f", data[44:48]))

if __name__ == "__main__":
    main()
