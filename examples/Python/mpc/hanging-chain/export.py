import os


def export_animation(ani):
    """Export the animation"""
    out = os.path.join(os.path.dirname(__file__), '..', '..', '..', '..',
                       'doxygen', 'sphinx', 'source', 'sphinxstatic',
                       'hanging-chain.html')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        f.write('<center>')
        f.write(ani.to_jshtml())
        f.write('</center>')
