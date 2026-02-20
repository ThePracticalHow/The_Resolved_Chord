#!/usr/bin/env python3
"""
LOTUS Documentation Validator
Checks for placeholder values and common issues in documentation.
"""

import os
import re
import sys
from pathlib import Path

def check_placeholders():
    """Check for TODO/FIXME placeholders in documentation."""
    issues = []

    # Files to check
    docs_files = [
        'README.md',
        'CONTRIBUTING.md',
        'CHANGELOG.md'
    ]

    placeholders = [
        r'TODO\.XXXXX',
        r'zenodo\.TODO',
        r'XXXX\.XXXXX',
        r'zenodo\.XXXXXXX'
    ]

    for doc_file in docs_files:
        if os.path.exists(doc_file):
            with open(doc_file, 'r', encoding='utf-8') as f:
                content = f.read()

            for placeholder in placeholders:
                if re.search(placeholder, content):
                    issues.append(f"Found placeholder '{placeholder}' in {doc_file}")

    return issues

def check_broken_links():
    """Check for obviously broken internal links."""
    issues = []

    if os.path.exists('README.md'):
        with open('README.md', 'r', encoding='utf-8') as f:
            content = f.read()

        # Check for markdown links to non-existent files
        link_pattern = r'\[([^\]]+)\]\(([^)]+)\)'
        links = re.findall(link_pattern, content)

        for text, url in links:
            # Skip external URLs
            if url.startswith(('http://', 'https://', 'mailto:')):
                continue

            # Check if local files exist
            if not os.path.exists(url):
                issues.append(f"Broken link: [{text}]({url})")

    return issues

def main():
    """Run all validation checks."""
    print("🔍 LOTUS Documentation Validator")
    print("=" * 40)

    all_issues = []

    # Check placeholders
    print("Checking for placeholders...")
    placeholder_issues = check_placeholders()
    all_issues.extend(placeholder_issues)

    # Check broken links
    print("Checking for broken links...")
    link_issues = check_broken_links()
    all_issues.extend(link_issues)

    # Report results
    if all_issues:
        print(f"\n❌ Found {len(all_issues)} issues:")
        for issue in all_issues:
            print(f"  - {issue}")
        print("\n💡 Fix these issues before publishing!")
        return 1
    else:
        print("\n✅ All checks passed!")
        return 0

if __name__ == '__main__':
    sys.exit(main())