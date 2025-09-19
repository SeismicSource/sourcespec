#!/usr/bin/env python3
"""
Script to reformat a BibTeX file:
- Add extra curly brackets to title only if not already double-braced
- Convert author first names to initials (adds dot if missing)
- Robust handling of nested braces, LaTeX escapes, and multi-line fields
"""
import shutil
import re
import argparse


def process_bibtex_entry(entry):
    """Process a single BibTeX entry string."""
    # --- Title processing ---
    def add_extra_curly(match):
        title_content = match.group(2).strip()
        comma = match.group(3)

        # Remove outer braces temporarily
        if title_content.startswith('{') and title_content.endswith('}'):
            inner = title_content[1:-1].strip()
            # Already double-braced?
            if inner.startswith('{') and inner.endswith('}'):
                return match.group(0)
            title_content = inner

        return f'{match.group(1)}{{{{{title_content}}}}}{comma}'

    entry = re.sub(
        r'(title\s*=\s*)({.*?})(,)',
        add_extra_curly,
        entry,
        flags=re.IGNORECASE | re.DOTALL
    )

    # --- Author processing ---
    def split_authors(authors_field):
        """
        Split authors by ' and ' ignoring 'and' inside braces.
        """
        authors = []
        brace_level = 0
        current = ''
        i = 0
        while i < len(authors_field):
            c = authors_field[i]
            if c == '{':
                brace_level += 1
            elif c == '}':
                brace_level -= 1
            # Only split on ' and ' at brace_level 0
            if (
                c == ' ' and
                authors_field[i:i+5] == ' and ' and
                brace_level == 0
            ):
                authors.append(current.strip())
                current = ''
                i += 4  # skip 'and '
            else:
                current += c
            i += 1
        if current.strip():
            authors.append(current.strip())
        return authors

    def initials(name):
        parts = name.split(',')
        if len(parts) == 2:
            last, firsts = [p.strip() for p in parts]
            new_firsts = []
            for n in firsts.split():
                if len(n) == 1:
                    new_firsts.append(f'{n}.')
                elif len(n) == 2 and n.endswith('.'):
                    new_firsts.append(n)
                else:
                    new_firsts.append(f'{n[0]}.')
            return f'{last}, {" ".join(new_firsts)}'
        return name

    # Regex captures:
    # group(1) = 'author = ', group(2) = {...}, group(3) = comma
    entry = re.sub(
        r'(author\s*=\s*){(.*?)}(,)',
        lambda m: f'{m.group(1)}{{' + ' and '.join(
            initials(n.strip()) for n in split_authors(m.group(2))
        ) + f'}}{m.group(3)}',
        entry,
        flags=re.IGNORECASE | re.DOTALL
    )

    return entry


def process_bibtex_file(input_file, output_file):
    """Process the entire BibTeX file."""
    with open(input_file, 'r', encoding='utf-8') as f:
        content = f.read()

    entries = re.split(r'@', content)[1:]
    processed_entries = [f'@{entry}' for entry in entries]
    processed_entries = [process_bibtex_entry(e) for e in processed_entries]

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(''.join(processed_entries))


def main():
    """Main function to handle argument parsing and file processing."""
    parser = argparse.ArgumentParser(description='Process a BibTeX file.')
    parser.add_argument(
        'input_file', help='Path to the input BibTeX file'
    )
    args = parser.parse_args()

    input_file = args.input_file
    backup_file = f'{input_file}.bak'

    # Create a backup
    shutil.copy2(input_file, backup_file)

    # Process and overwrite the original file
    process_bibtex_file(input_file, input_file)
    print(f'Original file backed up as {backup_file}')
    print(f'Processed BibTeX file updated in place: {input_file}')


if __name__ == '__main__':
    main()
