import sys


class Data:
    def __init__(self, fileName):
        print("Reading input file %s" % fileName)
        self.file = open(fileName, "r")
        self.loadData()
        self.file.close()

    def nextLine(self):
        line = self.file.readline().strip()
        if not line or line.startswith("#"):
            return self.nextLine()
        return line

    def loadData(self):
        self.vm_count = int(self.nextLine())  # VM的数量
        self.band_unitcost = float(self.nextLine())  # 带宽单位成本(η)

        self.vm_bandwidth = []  # VM需求的带宽
        for _ in range(self.vm_count):
            self.vm_bandwidth.append(int(self.nextLine()))

        self.vm_storage = []  # VM需求的存储
        for _ in range(self.vm_count):
            self.vm_storage.append(int(self.nextLine()))