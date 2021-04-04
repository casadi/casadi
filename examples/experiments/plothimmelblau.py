import matplotlib.pyplot as plt
import numpy as np 
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 13,
})
ext = '-2.png'

def himmelblau(x, y):
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

def meshplot(f, xlim, ylim, delta):
    x = np.arange(xlim[0], xlim[1], delta)
    y = np.arange(ylim[0], ylim[1], delta)
    X, Y = np.meshgrid(x, y)
    Z = f(X, Y)
    return X, Y, Z

xs = np.array([
3.3, 3,
2.9768120838964807, 2.356768229892948,
2.9825863590742134, 2.0714919737775643,
2.9725406720609122, 2.1000000000000001,
2.9724842368861419, 2.1000000000000001,
2.9724842353175864, 2.1000000000000001,
2.9724842353175864, 2.1000000000000001,
2.9724842353175864, 2.1000000000000001,
2.9724842353175864, 2.1000000000000001,
2.9724842353175864, 2.1000000000000001,
2.9724842353175864, 2.1000000000000001,
2.9724842353175864, 2.1000000000000001,
2.9724842353175864, 2.1000000000000001,
2.9724842353175864, 2.1000000000000001,
2.9724842353175864, 2.1000000000000001,
2.9724842353175864, 2.1000000000000001,
2.9724842353175864, 2.1000000000000001,
2.9724842353175864, 2.1000000000000001,
2.9724842353175864, 2.1000000000000001,
2.9724842353175864, 2.1000000000000001,
2.9724842353175864, 2.1000000000000001,
]).reshape(21, 2)
xs_gd = np.array([
3.3, 3,
3.05626, 2.6530999999999998,
2.9645602405429123, 2.4789242664181801,
2.9276011140376439, 2.3716567998060851,
2.9154439868141799, 2.298598584009429,
2.9152083498138426, 2.2454940935941861,
2.9205635546016735, 2.2050556958805769,
2.9282289676602042, 2.1731843649590354,
2.9364946766737052, 2.1474141830758828,
2.9444993458394761, 2.1261794247878587,
2.9518464548682286, 2.1084375172310512,
2.95839121654052, 2.1000000000000001,
2.9634594923706015, 2.1000000000000001,
2.9667131994272036, 2.1000000000000001,
2.968797184886542, 2.1000000000000001,
2.9701299888789632, 2.1000000000000001,
2.9709815668653734, 2.1000000000000001,
2.9715253401990958, 2.1000000000000001,
2.9718724302685353, 2.1000000000000001,
2.9720939224446825, 2.1000000000000001,
2.9722352430823382, 2.1000000000000001,
]).reshape(21, 2)
xs_lbfgs = np.array([
3.2999999999999998, 3,
3.05626, 2.6530999999999998,
2.8833557595352524, 2.3260813348768181,
2.8636362457906461, 2.1968743028552229,
2.8893863324619016, 2.1008521447169,
2.9360993482288813, 2.0633974971843023,
2.9958394724439454, 2.0986950752483673,
2.9729951800218388, 2.1019490101970502,
2.9718424792979659, 2.0998399050168235,
2.9724470226581992, 2.0999548614469652,
2.9724921002839166, 2.100004558430896,
2.9724845831850741, 2.1000002288193911,
2.9724842317048323, 2.0999999976263095,
2.9724842353169993, 2.0999999999999708,
2.9724842353175212, 2.0999999999999877,
2.9724842353175895, 2.0999999999999979,
2.9724842353175873, 2.0999999999999996,
2.9724842353175864, 2.1000000000000001,
2.9724842353175864, 2.1000000000000001,
2.9724842353175864, 2.1000000000000001,
2.9724842353175864, 2.1000000000000001,
]).reshape(21, 2)

plt.figure(figsize=(6.6, 8))
levels = np.arange(0, 54, 1)
X, Y, Z = meshplot(himmelblau, [2.5, 3.5], [1.9, 3.1], 0.025)
plt.contour(X, Y, Z, levels)
plt.plot(xs_gd[:,0], xs_gd[:,1], 'g.-', label='Gradient descent')
plt.plot(xs_lbfgs[:,0], xs_lbfgs[:,1], 'b.-', label='L-BFGS')
plt.plot(xs[:,0], xs[:,1], 'r.-', label='Newton')
plt.plot(3, 2, 'kx')
plt.axhline(2.1, linestyle=':', color='k', linewidth=1)
plt.legend()
plt.title('Himmelblau')
plt.gca().set_aspect('equal')
plt.savefig('himmelblau' + ext)

plt.figure(figsize=(6.6, 8))
levels = np.arange(0, 0.875, 0.0125)
xlim, ylim = [2.97224-0.05, 2.97224+0.05], [2.04, 2.16]
X, Y, Z = meshplot(himmelblau, xlim, ylim, 0.00025)
plt.contour(X, Y, Z, levels)
plt.plot(xs_gd[:,0], xs_gd[:,1], 'g.-', label='Gradient descent')
plt.plot(xs_lbfgs[:,0], xs_lbfgs[:,1], 'b.-', label='L-BFGS')
plt.plot(xs[:,0], xs[:,1], 'r.-', label='Newton')
plt.plot(3, 2, 'kx')
plt.axhline(2.1, linestyle=':', color='k', linewidth=1)
plt.legend()
plt.title('Himmelblau')
plt.xlim(xlim)
plt.ylim(ylim)
plt.gca().set_aspect('equal')
plt.savefig('himmelblau-close' + ext)

plt.figure(figsize=(6.6, 8))
h_end = himmelblau(*xs[-1])
plt.semilogy(abs(himmelblau(xs_gd[:,0], xs_gd[:,1]) - h_end), 'g.-', label='Gradient descent')
plt.semilogy(abs(himmelblau(xs_lbfgs[:,0], xs_lbfgs[:,1]) - h_end), 'b.-', label='L-BFGS')
plt.semilogy(abs(himmelblau(xs[:,0], xs[:,1]) - h_end), 'r.-', label='Newton')
plt.legend()
plt.title('Himmelblau convergence')
plt.ylabel(r'$\left| f(x) - f(x^\star) \right|$')
plt.savefig('himmelblau-conv-f' + ext)

plt.figure(figsize=(6.6, 8))
h_end = xs[-1]
plt.semilogy(np.linalg.norm(xs_gd - h_end, axis=-1), 'g.-', label='Gradient descent')
plt.semilogy(np.linalg.norm(xs_lbfgs - h_end, axis=-1), 'b.-', label='L-BFGS')
plt.semilogy(np.linalg.norm(xs - h_end, axis=-1), 'r.-', label='Newton')
plt.legend()
plt.title('Himmelblau convergence')
plt.ylabel(r'$\left\|x - x^\star \right\|_2$')
plt.savefig('himmelblau-conv-x' + ext)

plt.show()