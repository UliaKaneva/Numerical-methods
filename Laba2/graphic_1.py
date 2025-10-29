import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(-2, 2, 400)

y = np.sin(x) - 2 * x ** 2 + 0.5
x_1 = np.linspace(0.7, 0.85, 100)
phi = ((np.sin(x_1) + 0.5) / 2) ** 0.5
phi_pr = np.cos(x_1) / (2 * (2 * np.sin(x_1) - 1) ** 0.5)

plt.figure(figsize=(10, 6))

plt.plot(x, y, label='$y = sin(x) - 2x^2 + 0.5$', color='red', linewidth=2)
plt.plot(x_1, phi, label='$ \\varphi(x) = \\sqrt{\\frac{sin(x) + 0.5}{2}}$', color='blue', linewidth=2)
plt.plot(x_1, phi_pr, label='$\\varphi\'(x) = \\frac{cos(x_1)}{2 \cdot \\sqrt{2 \cdot sin(x_1) - 1}}$', color='black', linewidth=2)

plt.title('График функции $y = sin(x) - 2 * x^2 + 0.5$')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True, alpha=0.3)
plt.axhline(y=0, color='k', linewidth=0.5)
plt.axvline(x=0, color='k', linewidth=0.5)

plt.show()
