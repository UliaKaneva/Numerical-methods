import matplotlib.pyplot as plt
import numpy as np

x_2 = np.linspace(-1, 3, 400)

y_1 = 1 + np.cos(x_2)
x_1 = np.linspace(-1, 3, 400)
y_2 = 1 + np.sin(x_1)

plt.figure(figsize=(10, 6))

plt.plot(x_1, y_2, label='$x_2 = 1 + sin(x_1)$', color='red', linewidth=2)
plt.plot(y_1, x_2, label='$x_1 = 1 + cos(x_2)$', color='blue', linewidth=2)

plt.title('График системы')
plt.xlabel('$x_1$')
plt.ylabel('$x_2$')
plt.legend()
plt.grid(True, alpha=0.3)
plt.axhline(y=0, color='k', linewidth=0.5)
plt.axvline(x=0, color='k', linewidth=0.5)

plt.show()
