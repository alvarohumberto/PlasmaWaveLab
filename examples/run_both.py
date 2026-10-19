"""Run both reproducible command-line examples from an installed checkout."""
from pathlib import Path
from plasmawavelab.__main__ import main
import sys

if __name__ == '__main__':
    for model in ('sagdeev','kdv'):
        sys.argv = ['plasmawavelab', model, '--output', str(Path('results') / model)]
        main()
