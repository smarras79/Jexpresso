"""
PDF Parser for extracting equation information from PDF files.
"""
import pdfplumber
from typing import Dict, List, Optional
from pathlib import Path


class PDFParser:
    """Extracts text and images from PDF files."""

    def __init__(self, pdf_path: str):
        """
        Initialize PDF parser.

        Args:
            pdf_path: Path to the PDF file
        """
        self.pdf_path = Path(pdf_path)
        if not self.pdf_path.exists():
            raise FileNotFoundError(f"PDF file not found: {pdf_path}")

    def extract_text(self) -> str:
        """
        Extract all text from the PDF.

        Returns:
            Concatenated text from all pages
        """
        text_content = []

        with pdfplumber.open(self.pdf_path) as pdf:
            for page_num, page in enumerate(pdf.pages, 1):
                text = page.extract_text()
                if text:
                    text_content.append(f"--- Page {page_num} ---\n{text}")

        return "\n\n".join(text_content)

    def extract_images(self, output_dir: Optional[str] = None) -> List[Dict]:
        """
        Extract images from the PDF.

        Args:
            output_dir: Directory to save extracted images (optional)

        Returns:
            List of image metadata dictionaries
        """
        images = []

        with pdfplumber.open(self.pdf_path) as pdf:
            for page_num, page in enumerate(pdf.pages, 1):
                page_images = page.images
                for img_num, img in enumerate(page_images):
                    images.append({
                        'page': page_num,
                        'image_num': img_num,
                        'bbox': (img['x0'], img['top'], img['x1'], img['bottom']),
                        'width': img['width'],
                        'height': img['height']
                    })

        return images

    def get_page_count(self) -> int:
        """Get the number of pages in the PDF."""
        with pdfplumber.open(self.pdf_path) as pdf:
            return len(pdf.pages)
