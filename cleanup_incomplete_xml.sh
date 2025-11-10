#!/bin/bash
# Script to remove incomplete BLAST XML files (those missing </BlastOutput> closing tag)

echo "Checking for incomplete XML files..."

# Find and count incomplete files
incomplete_count=$(find outputs/*/phase2_synteny -name "*.xml" -exec sh -c 'tail -1 "$1" 2>/dev/null | grep -q "</BlastOutput>" || echo "$1"' sh {} \; 2>/dev/null | wc -l)

echo "Found $incomplete_count incomplete XML files"

if [ $incomplete_count -gt 0 ]; then
    echo "Creating backup list of files to be deleted..."
    find outputs/*/phase2_synteny -name "*.xml" -exec sh -c 'tail -1 "$1" 2>/dev/null | grep -q "</BlastOutput>" || echo "$1"' sh {} \; 2>/dev/null > incomplete_xml_files.txt

    echo "First 10 files to be deleted:"
    head -10 incomplete_xml_files.txt

    echo ""
    echo "Total files to delete: $(wc -l < incomplete_xml_files.txt)"
    echo "Press Enter to delete these files, or Ctrl+C to cancel..."
    read

    # Delete the incomplete files
    while IFS= read -r file; do
        rm -f "$file"
    done < incomplete_xml_files.txt

    echo "Deleted $incomplete_count incomplete XML files"
    echo "List saved to incomplete_xml_files.txt for reference"
else
    echo "No incomplete XML files found!"
fi

# Check what remains
remaining=$(find outputs/*/phase2_synteny -name "*.xml" | wc -l)
echo "Remaining XML files: $remaining"