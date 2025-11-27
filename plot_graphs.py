import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv
import sys

def read_data(filename):
    data = {'MPI': {}, 'OpenMP': {}, 'Hybrid2': {}, 'Hybrid4': {}}
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            next(reader, None)

            for row in reader:
                if not row or len(row) < 3: continue
                dtype = row[0].strip()
                try:
                    val = int(row[1].strip())
                    time = float(row[2].strip())
                    if dtype in data:
                        data[dtype][val] = time
                except ValueError: continue
    except FileNotFoundError:
        print(f"Файл {filename} не найден.")
        sys.exit(1)
    return data

def get_metrics(dataset, base_time):
    if not dataset: return [], [], [], []
    x = sorted(dataset.keys())
    t = [dataset[i] for i in x]
    # Speedup = BaseTime / CurrentTime
    speedup = [base_time / val for val in t]
    # Efficiency = Speedup / (ImprovementFactor)
    efficiency = [s / c for s, c in zip(speedup, x)]
    return x, t, speedup, efficiency

def plot_graphs(input_file, output_file, title_suffix):
    d = read_data(input_file)
    if not any(d.values()): return
    t_base = d['MPI'].get(1, list(d['OpenMP'].values())[0] if d['OpenMP'] else 1.0)
    mx, mt, ms, me = get_metrics(d['MPI'], t_base)
    ox, ot, os, oe = get_metrics(d['OpenMP'], t_base)

    # Hybrid2 (X=threads): Speedup считаем относительно MPI=2
    base_h2 = d['MPI'].get(2, t_base)
    h2x, h2t, h2s, h2e = get_metrics(d['Hybrid2'], base_h2)

    # Hybrid4 (X=threads): Speedup считаем относительно MPI=4
    base_h4 = d['MPI'].get(4, t_base)
    h4x, h4t, h4s, h4e = get_metrics(d['Hybrid4'], base_h4)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f'Результаты: {title_suffix}', fontsize=14)

    styles = [
        ('MPI (Processes)', mx, mt, ms, me, 'o-', '#1f77b4'),
        ('OpenMP (Threads)', ox, ot, os, oe, 's-', '#d62728'),
        ('Hybrid (2 Proc x N Threads)', h2x, h2t, h2s, h2e, '^-', '#2ca02c'),
        ('Hybrid (4 Proc x N Threads)', h4x, h4t, h4s, h4e, 'D-', '#9467bd')
    ]

    # 1. Время
    for label, x, t, s, e, fmt, clr in styles:
        if x: ax1.plot(x, t, fmt, label=label, color=clr, linewidth=2)

    ax1.set_title('Время выполнения')
    ax1.set_xlabel('N)')
    ax1.set_ylabel('Время (с)')
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend()

    # 2. Ускорение
    for label, x, t, s, e, fmt, clr in styles:
        if x:
            l_s = label
            if "Hybrid" in label: l_s += " (vs pure MPI)"
            ax2.plot(x, s, fmt, label=l_s, color=clr, linewidth=2)

    ax2.set_title('Ускорение')
    ax2.set_xlabel('N')
    ax2.set_ylabel('Кратность')
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.legend()

    # 3. Эффективность
    for label, x, t, s, e, fmt, clr in styles:
        if x: ax3.plot(x, e, fmt, label=label, color=clr, linewidth=2)

    ax3.set_title('Эффективность')
    ax3.set_xlabel('N')
    ax3.set_ylabel('Коэф.')
    ax3.set_ylim(0, 1.2)
    ax3.grid(True, linestyle='--', alpha=0.6)
    # ax3.legend()

    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"Сохранен: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 4: sys.exit(1)
    plot_graphs(sys.argv[1], sys.argv[2], sys.argv[3])