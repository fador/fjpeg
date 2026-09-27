import numpy as np
import os

def create_images():
    w, h = 640, 480
    yy, xx = np.mgrid[0:h, 0:w].astype(np.float32)

    # 1. Gradient + Smooth Texture
    y1 = np.clip((xx / w * 180 + yy / h * 50) + 20 * np.sin(xx / 10.0) * np.cos(yy / 10.0) + 15 * np.sin(xx * yy / 500.0), 16, 235).astype(np.uint8)
    cb1 = np.clip(128 + 30 * np.sin(xx[::2, ::2] / 20.0), 16, 240).astype(np.uint8)
    cr1 = np.clip(128 + 30 * np.cos(yy[::2, ::2] / 20.0), 16, 240).astype(np.uint8)
    with open('test_grad_texture.yuv', 'wb') as f:
        f.write(y1.tobytes())
        f.write(cb1.tobytes())
        f.write(cr1.tobytes())

    # 2. Natural-like photographic scene (landscape with sky, mountains, grass texture)
    sky = 180.0 - 50.0 * (yy / (h * 0.45))
    sky = np.clip(sky, 130, 220)
    # Mountains
    mountain_h = 0.45 * h + 40.0 * np.sin(xx / 45.0) + 25.0 * np.sin(xx / 18.0) + 10.0 * np.cos(xx / 7.0)
    mountain_mask = (yy >= mountain_h)
    mountain_tex = 70.0 + 30.0 * np.sin((xx + yy) / 12.0) + 15.0 * np.cos((xx * 2 - yy) / 15.0)
    # Foreground grass / trees
    fg_h = 0.70 * h + 20.0 * np.sin(xx / 30.0)
    fg_mask = (yy >= fg_h)
    np.random.seed(42)
    grass_tex = 50.0 + 25.0 * np.sin(xx / 3.0) * np.cos(yy / 3.0) + np.random.normal(0, 8.0, (h, w))
    
    y2 = np.where(fg_mask, grass_tex, np.where(mountain_mask, mountain_tex, sky))
    y2 = np.clip(y2, 16, 235).astype(np.uint8)
    
    cb2 = np.where(fg_mask[::2, ::2], 110, np.where(mountain_mask[::2, ::2], 125, 140)).astype(np.uint8)
    cr2 = np.where(fg_mask[::2, ::2], 115, np.where(mountain_mask[::2, ::2], 130, 120)).astype(np.uint8)
    with open('test_natural.yuv', 'wb') as f:
        f.write(y2.tobytes())
        f.write(cb2.tobytes())
        f.write(cr2.tobytes())

    # 3. Geometric / sharp edges
    y3 = np.zeros((h, w), dtype=np.float32) + 200.0
    for i in range(0, w, 40):
        y3[:, i:i+20] -= 120.0
    for j in range(0, h, 40):
        y3[j:j+20, :] -= 40.0
    y3 = np.clip(y3, 16, 235).astype(np.uint8)
    cb3 = np.full((h//2, w//2), 128, dtype=np.uint8)
    cr3 = np.full((h//2, w//2), 128, dtype=np.uint8)
    with open('test_edges.yuv', 'wb') as f:
        f.write(y3.tobytes())
        f.write(cb3.tobytes())
        f.write(cr3.tobytes())

    print('Generated test_grad_texture.yuv, test_natural.yuv, test_edges.yuv')

if __name__ == '__main__':
    create_images()
