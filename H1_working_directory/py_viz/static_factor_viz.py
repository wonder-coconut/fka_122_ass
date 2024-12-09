import matplotlib.pyplot as plt
import numpy as np

#sq output read
sq = open("../op_text/sq_t7.txt").read().split('\n')
sq.pop()
sq = np.array([float(item) for item in sq])

#wavevector grid read
q_grid = open("../op_text/qgrid_t7.txt").read().split('\n')
q_grid.pop()
for i in range(len(q_grid)):
    q_grid[i] = q_grid[i].split("\t")
    q_grid[i] = [float(item) for item in q_grid[i]]

q_grid = np.array(q_grid)

q_res = np.array([np.linalg.norm(q) for q in q_grid])
q_0_ind = np.where(q_res == 0)

q_res = np.delete(q_res, q_0_ind)
sq = np.delete(sq, q_0_ind)

plt.plot(sq)
plt.show()