import subprocess, os, sys
import numpy as np

def bdrint(rate, dist):
    p = np.polyfit(dist, np.log(rate), 3)
    return np.polyint(p)

def bd_rate(r1, d1, r2, d2):
    min_d = max(min(d1), min(d2))
    max_d = min(max(d1), max(d2))
    if min_d >= max_d:
        return 0.0
    p1 = bdrint(r1, d1)
    p2 = bdrint(r2, d2)
    int1 = np.polyval(p1, max_d) - np.polyval(p1, min_d)
    int2 = np.polyval(p2, max_d) - np.polyval(p2, min_d)
    avg_diff = (int2 - int1) / (max_d - min_d)
    return (np.exp(avg_diff) - 1.0) * 100.0

def psnr(orig, dec):
    mse = np.mean((orig.astype(np.float64) - dec.astype(np.float64))**2)
    if mse == 0: return 99.99
    return 10.0 * np.log10(255.0**2 / mse)

def evaluate_image(yuv_file, width=640, height=480, fjpeg_bin='build/Release/fjpeg.exe'):
    with open(yuv_file, 'rb') as f:
        orig_y = np.frombuffer(f.read(width * height), dtype=np.uint8)
    
    qualities = [30, 50, 70, 85]
    configs = [
        ('baseline', []),
        ('trellis', ['-t']),
        ('progressive', ['-p']),
        ('prog+trellis', ['-p', '-t']),
        ('arithmetic', ['-a']),
        ('arith+trellis', ['-a', '-t']),
    ]
    
    results = {}
    for name, flags in configs:
        rates = []
        psnrs = []
        for q in qualities:
            jpg = f'tmp_{name}_{q}.jpg'
            yuv_dec = f'tmp_{name}_{q}.yuv'
            cmd = [fjpeg_bin, '-i', yuv_file, '-r', f'{width}x{height}', '-q', str(q), '-o', jpg] + flags
            ret = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if ret.returncode != 0:
                print(f"Error encoding {jpg}")
                continue
            sz = os.path.getsize(jpg)
            
            # decode with fjpeg
            cmd_dec = [fjpeg_bin, '-d', '-i', jpg, '-o', yuv_dec]
            subprocess.run(cmd_dec, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            with open(yuv_dec, 'rb') as f:
                dec_y = np.frombuffer(f.read(width * height), dtype=np.uint8)
            p = psnr(orig_y, dec_y)
            rates.append(sz)
            psnrs.append(p)
            try:
                os.remove(jpg)
                os.remove(yuv_dec)
            except:
                pass
        results[name] = (rates, psnrs)
    
    base_rates, base_psnrs = results['baseline']
    print(f"\nResults for {yuv_file}:")
    print(f"{'Config':<15} | {'q=30 (B / dB)':<18} | {'q=50 (B / dB)':<18} | {'q=70 (B / dB)':<18} | {'q=85 (B / dB)':<18} | {'BD-rate vs Base':<16}")
    print("-" * 115)
    for name, (rates, psnrs) in results.items():
        q_strs = [f"{rates[i]} / {psnrs[i]:.2f}" for i in range(4)]
        if name == 'baseline':
            bd_str = "0.00% (ref)"
        else:
            bd = bd_rate(base_rates, base_psnrs, rates, psnrs)
            bd_str = f"{bd:+.2f}%"
        print(f"{name:<15} | {q_strs[0]:<18} | {q_strs[1]:<18} | {q_strs[2]:<18} | {q_strs[3]:<18} | {bd_str:<16}")

    # Lossless benchmark
    jpg_ll = 'tmp_lossless.jpg'
    yuv_ll = 'tmp_lossless.yuv'
    subprocess.run([fjpeg_bin, '-i', yuv_file, '-r', f'{width}x{height}', '-lossless', '-o', jpg_ll], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    ll_sz = os.path.getsize(jpg_ll)
    subprocess.run([fjpeg_bin, '-d', '-i', jpg_ll, '-o', yuv_ll], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open(yuv_file, 'rb') as f1, open(yuv_ll, 'rb') as f2:
        exact = (f1.read() == f2.read())
    raw_sz = width * height * 3 // 2
    ratio = raw_sz / ll_sz
    bpp = (ll_sz * 8) / (width * height)
    print(f"Lossless mode   | Compressed: {ll_sz} bytes | Ratio: {ratio:.2f}:1 | {bpp:.2f} bpp | Exact bit-for-bit: {exact}")
    try:
        os.remove(jpg_ll)
        os.remove(yuv_ll)
    except:
        pass

    return results

def evaluate_12bit_image(yuv_file, width=640, height=480, fjpeg_bin='build/Release/fjpeg.exe'):
    if not os.path.exists(fjpeg_bin):
        fjpeg_bin = 'build/Debug/fjpeg.exe'
    print(f"\n--- 12-bit Benchmark for {yuv_file} ---")
    raw_sz = width * height * 3 // 2 * 2  # 16-bit words

    # 1. 12-bit Lossless
    jpg_ll = 'tmp_12bit_lossless.jpg'
    yuv_ll = 'tmp_12bit_lossless.yuv'
    subprocess.run([fjpeg_bin, '-i', yuv_file, '-r', f'{width}x{height}', '-b', '12', '-lossless', '-o', jpg_ll], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    ll_sz = os.path.getsize(jpg_ll)
    subprocess.run([fjpeg_bin, '-d', '-i', jpg_ll, '-o', yuv_ll], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open(yuv_file, 'rb') as f1, open(yuv_ll, 'rb') as f2:
        exact = (f1.read() == f2.read())
    ratio = raw_sz / ll_sz
    bpp = (ll_sz * 8) / (width * height)
    print(f"12-bit Lossless  | Compressed: {ll_sz} bytes | Ratio: {ratio:.2f}:1 | {bpp:.2f} bpp | Exact bit-for-bit: {exact}")

    # 2. 12-bit DCT
    with open(yuv_file, 'rb') as f:
        orig = np.frombuffer(f.read(), dtype=np.uint16)
    for q in [50, 70, 90]:
        jpg_dct = f'tmp_12bit_q{q}.jpg'
        yuv_dct = f'tmp_12bit_q{q}.yuv'
        jpg_arith = f'tmp_12bit_arith_q{q}.jpg'
        subprocess.run([fjpeg_bin, '-i', yuv_file, '-r', f'{width}x{height}', '-b', '12', '-q', str(q), '-o', jpg_dct], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([fjpeg_bin, '-i', yuv_file, '-r', f'{width}x{height}', '-b', '12', '-a', '-q', str(q), '-o', jpg_arith], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        dct_sz = os.path.getsize(jpg_dct)
        arith_sz = os.path.getsize(jpg_arith)
        subprocess.run([fjpeg_bin, '-d', '-i', jpg_dct, '-o', yuv_dct], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        with open(yuv_dct, 'rb') as f:
            dec = np.frombuffer(f.read(), dtype=np.uint16)
        mse = np.mean((orig.astype(float) - dec.astype(float))**2)
        p = 10.0 * np.log10(4095.0**2 / mse) if mse > 0 else 99.99
        saving = (dct_sz - arith_sz) / dct_sz * 100.0
        print(f"12-bit DCT q={q:<2}  | Base: {dct_sz:6d} B | Arith: {arith_sz:6d} B ({saving:+.1f}%) | PSNR: {p:.2f} dB")
        try:
            os.remove(jpg_dct)
            os.remove(jpg_arith)
            os.remove(yuv_dct)
        except:
            pass

    try:
        os.remove(jpg_ll)
        os.remove(yuv_ll)
    except:
        pass

if __name__ == '__main__':
    bin_path = 'build/Release/fjpeg.exe' if os.path.exists('build/Release/fjpeg.exe') else 'build/Debug/fjpeg.exe'
    for img in ['test_grad_texture.yuv', 'test_natural.yuv', 'test_edges.yuv']:
        if not os.path.exists(img):
            import generate_test_images
            generate_test_images.create_images()
        evaluate_image(img, fjpeg_bin=bin_path)

    if os.path.exists('test_natural_12bit.yuv'):
        evaluate_12bit_image('test_natural_12bit.yuv', fjpeg_bin=bin_path)

