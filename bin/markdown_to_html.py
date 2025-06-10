#!/usr/bin/env python3

import argparse
import sys
import markdown

def convert_markdown(md_file):
    """Convert markdown to HTML with basic extensions only."""
    with open(md_file, 'r') as f:
        md_content = f.read()

    # Use only basic extensions that come with markdown
    html = markdown.markdown(
        md_content,
        extensions=[
            'markdown.extensions.tables',
            'markdown.extensions.fenced_code',
            'markdown.extensions.codehilite',
            'markdown.extensions.toc'
        ]
    )
    return html

def wrap_html(body_html, title="Analysis Results"):
    """Wrap the converted markdown in a complete HTML document."""
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            line-height: 1.6;
            color: #333;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
        }}
        h1, h2, h3 {{ color: #2c3e50; }}
        code {{
            background-color: #f4f4f4;
            padding: 2px 4px;
            border-radius: 3px;
            font-family: 'Courier New', monospace;
        }}
        pre {{
            background-color: #f4f4f4;
            padding: 10px;
            border-radius: 5px;
            overflow-x: auto;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            margin: 20px 0;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 8px;
            text-align: left;
        }}
        th {{ background-color: #f2f2f2; }}
        .toc {{ background-color: #f9f9f9; padding: 15px; border-radius: 5px; }}
    </style>
</head>
<body>
{body_html}
</body>
</html>"""

def main():
    parser = argparse.ArgumentParser(description='Convert Markdown to HTML')
    parser.add_argument('mdfile', type=argparse.FileType('r'), help='Input markdown file')
    parser.add_argument('-o', '--output', required=True, help='Output HTML file')
    parser.add_argument('--title', default='Analysis Results', help='HTML page title')

    args = parser.parse_args()

    try:
        converted_md = convert_markdown(args.mdfile.name)
        full_html = wrap_html(converted_md, args.title)

        with open(args.output, 'w') as f:
            f.write(full_html)

        print(f"Successfully converted {args.mdfile.name} to {args.output}")

    except Exception as e:
        print(f"Error converting markdown: {e}", file=sys.stderr)
        return 1

    finally:
        args.mdfile.close()

    return 0

if __name__ == '__main__':
    sys.exit(main())

