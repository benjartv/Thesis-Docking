torsion = []
idatm = []
atm = []
x = []
y = []
z = []
occupancy = []
temp = []
charge = []
atmType = []


def readBranch(content, a, atm1, atm2):
	i = a
	branch = []
	while "ENDBRANCH" != content[i][0]:
		if content[i][0] == "HETATM":
			branch.append(int(content[i][1]))
			
			idatm.append(content[i][1])
			atm.append(content[i][2])
			x.append(content[i][6])
			y.append(content[i][7])
			z.append(content[i][8])
			occupancy.append(content[i][9])
			temp.append(content[i][10])
			charge.append(content[i][11])
			atmType.append(content[i][12])
			
			i+=1
		elif content[i][0] == "BRANCH":
			branchdata = readBranch(content, i+1, content[1], content[2])
			branch.append((int(content[i][1]), int(content[i][2])))
			branch.append(branchdata[0])
			i = branchdata[1]
	return [branch, i]



file = open("lig.pdbqt", "r")


data = file.readlines()
file.close()
content = []
for line in data:
	content.append(line.strip().split())

root = []
i = 1
while True:
	if content[i][0] == "ENDROOT":
		break
	idatm.append(content[i][1])
	atm.append(content[i][2])
	x.append(content[i][6])
	y.append(content[i][7])
	z.append(content[i][8])
	occupancy.append(content[i][9])
	temp.append(content[i][10])
	charge.append(content[i][11])
	atmType.append(content[i][12])
	root.append(content[i][1])
	i+=1



branch = []
while content[i][0] != "TORSDOF":
	if content[i][0] == "BRANCH":
		branch.append([content[i][1],content[i][2]])
	elif content[i][0] == "HETATM":
		idatm.append(content[i][1])
		atm.append(content[i][2])
		x.append(content[i][6])
		y.append(content[i][7])
		z.append(content[i][8])
		occupancy.append(content[i][9])
		temp.append(content[i][10])
		charge.append(content[i][11])
		atmType.append(content[i][12])
	i+=1

for link in branch:
	print link

print idatm
print atm



'''
branch = []
while content[i][0] != "TORSDOF":
	if content[i][0] == "BRANCH":
		newbranchdata = readBranch(content, i+1, content[i][1], content[i][2])
		branch.append((int(content[i][1]), int(content[i][2])))
		branch.append(newbranchdata[0])
		i = newbranchdata[1]
	else:
		i+=1

print branch




root = []
while True:
	if "ENDROOT" in content:
		break
	idatm.append(content[1])
	atm.append(content[2])
	x.append(content[6])
	y.append(content[7])
	z.append(content[8])
	occupancy.append(content[9])
	temp.append(content[10])
	charge.append(content[11])
	atmType.append(content[12])
	root.append(content[1])
idb = -1
branch = []
while "TORSDOF" not in content:
	content = file.readline().strip().split()
	if content[0] == "BRANCH":
		batm1 = content[1]
		batm2 = content[2]
		branch.append(readBranch(file, batm1,batm2))
		
		idb += 1
		branch[idb].append((content[1],content[2]))
	elif content[0] == "HETATM":
		idatm.append(content[1])
		atm.append(content[2])
		x.append(content[6])
		y.append(content[7])
		z.append(content[8])
		occupancy.append(content[9])
		temp.append(content[10])
		charge.append(content[11])
		atmType.append(content[12])
		branch[idb].append(content[1])
		

file.close()
for atm in root:
	print atm

print branch
'''