import math
from math import fabs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def Boundary_condition(a, b) :
	if a == 0.0 * dx :
		return 1
	elif a <= 57.0 * dx :
		O1 = 0.0 * dx
		O2 = O1 + 160 * dx
		O3 = O2 + 160 * dx
		O4 = O3 + 160 * dx
		width = (66 - (a / dx + int(a / dx / 6) - 1)) * dx
		if (width / 2.0 - fabs(b - O1)) > -1e-6 or (width / 2.0 - fabs(b - O2)) > -1e-6 or (width / 2.0 - fabs(b - O3)) > -1e-6 or (width / 2.0 - fabs(b - O4)) > -1e-6 :
			return 1
		else :
			return 0
	elif a >= 78 * dx and a <= 98 * dx :
		O1 = 80.0 * dx
		O2 = O1 + 160.0 * dx
		O3 = O2 + 160.0 * dx
		O4 = O3 + 160.0 * dx
		width = 70.0 * dx
		if (width / 2.0 - fabs(b - O1)) > -1e-6 or (width / 2.0 - fabs(b - O2)) > -1e-6 or (width / 2.0 - fabs(b - O3)) > -1e-6 or (width / 2.0 - fabs(b - O4)) > -1e-6 :
			return 2
		else :
			return 0
	else :
		return 0

def Potential() :
	global phi
	for j in range(Cell_y):
		for i in range(Cell_x):
			if Boundary_condition(i * dx, j * dx) == 0:
				phi[i][j] = 0.0
			elif Boundary_condition(i * dx, j * dx) == 1:
				phi[i][j] = P
			elif  Boundary_condition(i * dx, j * dx) == 2:
				phi[i][j] = 0.0
	phi_pre = phi.copy()
	dphi_max = 100
	pi = math.pi
	aL = aR = 1 / (dx ** 2)
	aU = aD = 1 / (dx ** 2)
	aP = 4 / (dx ** 2)
	ω = 2 / (1 + (1 - ((math.cos(pi / Cell_x) + math.cos(pi / Cell_y)) / 2) ** 2) ** 0.5)
	while dphi_max > 0.1:
		dphi_max = 0
		for j in range(Cell_y):
			for i in range(Cell_x):
				if Boundary_condition(i * dx, j * dx) != 0:
					continue
				elif i == Cell_x - 1 and j == 0:
					phi[i][j] = (1 - ω) * phi_pre[i][j] + ω * ((2 * aL * phi[i - 1][j] + 2 * aU * phi_pre[i][j + 1]) / aP + e / ep0 * N_ca[i][j] / aP)
				elif i == Cell_x - 1 and j == Cell_y - 1:
					phi[i][j] = (1 - ω) * phi_pre[i][j] + ω * ((2 * aL * phi[i - 1][j] + 2 * aD * phi[i][j - 1]) / aP + e / ep0 * N_ca[i][j] / aP)
				elif j == 0:
					phi[i][j] = (1 - ω) * phi_pre[i][j] + ω * ((aL * phi[i - 1][j] + aR * phi_pre[i + 1][j] + 2 * aU * phi_pre[i][j + 1]) / aP + e / ep0 * N_ca[i][j] / aP)
				elif j == Cell_y - 1:
					phi[i][j] = (1 - ω) * phi_pre[i][j] + ω * ((aL * phi[i - 1][j] + aR * phi_pre[i + 1][j] + 2 * aD * phi[i][j - 1]) / aP + e / ep0 * N_ca[i][j] / aP)
				elif i == Cell_x - 1:
					phi[i][j] = (1 - ω) * phi_pre[i][j] + ω * ((2 * aL * phi[i - 1][j] + aD * phi[i][j - 1] + aU * phi_pre[i][j + 1]) / aP + e / ep0 * N_ca[i][j] / aP)
				else:
					phi[i][j] = (1 - ω) * phi_pre[i][j] + ω * ((aL * phi[i - 1][j] + aR * phi_pre[i + 1][j] + aD * phi[i][j - 1] + aU * phi_pre[i][j + 1]) / aP + e / ep0 * N_ca[i][j] / aP)
				dphi = fabs(phi[i][j] - phi_pre[i][j])
				if dphi > dphi_max:
					dphi_max = dphi
		phi_pre = phi.copy()
		print("dphi_max=", dphi_max)

def Plot():
	x = np.linspace(0, Cell_x * dx, Cell_x)
	y = np.linspace(0, Cell_y * dx, Cell_y)
	Y, X = np.meshgrid(y, x)																#编织网格
	levels = MaxNLocator(nbins = P / 100).tick_values(0, P)									#设置参照表
	plt.contourf(X, Y, phi, levels = levels, cmap = plt.cm.jet)
	cbar = plt.colorbar()
	cbar.set_label('Unit: V', rotation = -90,va = 'bottom', fontsize = 20)							#设置标题
	cbar.ax.tick_params(labelsize = 10)														#设置轴上字体大小
	plt.show()


Cell_x = 601
Cell_y = 561
dx = 5e-6
e = 1.6e-19
ep0 = 8.854e-12
P = 2100.0

phi = np.array([[0.0] * Cell_y for i in range(Cell_x)])
N_ca = phi.copy()

Potential()
Plot()
