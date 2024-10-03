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
        self.vm_count = int(self.nextLine())  # VM的数量(0~n-1)
        self.pm_count = int(self.nextLine())  # PM的数量(1~n)
        self.cluster_count = int(self.nextLine())  # 集群的数量(1~n)
        self.band_unitcost = float(self.nextLine())  # 带宽单位成本(符号:η)
        self.elasticity_lb = float(self.nextLine())  # 负载弹性下限(符号:τ)

        self.vm_mycluster = []  # VM所属的集群
        for _ in range(self.vm_count):
            self.vm_mycluster.append(int(self.nextLine()))

        self.vm_mypm = []  # VM所属的物理机
        for _ in range(self.vm_count):
            self.vm_mypm.append(int(self.nextLine()))
        
        self.pm_CPUcore = []  # 物理机CPU核心数o(mi)
        for _ in range(self.pm_count):
            self.pm_CPUcore.append(int(self.nextLine()))

        self.pm_storage = []  # 物理机存储容量r(mi)
        for _ in range(self.pm_count):
            self.pm_storage.append(float(self.nextLine()))

        # 注意此处的带宽指公网带宽，仅作为一项资源，集群内vm间的通信已由“数据量×跳数×单位开销”表示
        self.pm_bw = []  # 物理机带宽b(mi)
        for _ in range(self.pm_count):
            self.pm_bw.append(float(self.nextLine()))

        self.vm_CPUuti = []  # VM需求的CPU利用率c(vk)
        for _ in range(self.vm_count):
            self.vm_CPUuti.append(float(self.nextLine()) / 100)

        self.vm_CPUcore = []  # VM需求的CPU核心o(vk)
        for _ in range(self.vm_count):
            self.vm_CPUcore.append(int(self.nextLine()))

        self.vm_storage = []  # VM需求的存储r(vk)
        for _ in range(self.vm_count):
            self.vm_storage.append(float(self.nextLine()))

        # 注意事项同上
        self.vm_bw = []  # VM需求的带宽b(vk)
        for _ in range(self.vm_count):
            self.vm_bw.append(float(self.nextLine()))

        self.hop_matrix = {}  # 跳数矩阵
        for i in range(self.pm_count):
            line = self.nextLine()
            for j, value in enumerate(line.split()):
                self.hop_matrix[(i, j)] = int(value)

        self.cluster_innerpakage_size = []  # 集群内数据包大小
        for _ in range(self.cluster_count):
            self.cluster_innerpakage_size.append(float(self.nextLine()))
