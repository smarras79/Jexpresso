"""
PDF Parser for extracting text from PDF files.
"""

using PDFIO

"""
    parse_pdf(pdf_path::String) -> String

Extract all text from a PDF file.

# Arguments
- `pdf_path::String`: Path to the PDF file

# Returns
- Concatenated text from all pages

# Example
```julia
text = parse_pdf("equations.pdf")
```
"""
function parse_pdf(pdf_path::String)
    if !isfile(pdf_path)
        error("PDF file not found: $pdf_path")
    end

    text_content = String[]

    try
        doc = pdDocOpen(pdf_path)
        npages = pdDocGetPageCount(doc)

        for page_num in 1:npages
            page = pdDocGetPage(doc, page_num)
            text = pdPageExtractText(IOBuffer(), page)
            page_text = String(take!(text))

            if !isempty(strip(page_text))
                push!(text_content, "--- Page $page_num ---\n$page_text")
            end
        end

        pdDocClose(doc)
    catch e
        @warn "Error parsing PDF with PDFIO: $e"
        @info "Falling back to simple text extraction..."

        # Fallback: try to read as text
        try
            text = read(pdf_path, String)
            push!(text_content, text)
        catch
            error("Could not extract text from PDF: $pdf_path")
        end
    end

    return join(text_content, "\n\n")
end

"""
    get_page_count(pdf_path::String) -> Int

Get the number of pages in a PDF file.

# Arguments
- `pdf_path::String`: Path to the PDF file

# Returns
- Number of pages

# Example
```julia
npages = get_page_count("equations.pdf")
```
"""
function get_page_count(pdf_path::String)
    if !isfile(pdf_path)
        error("PDF file not found: $pdf_path")
    end

    try
        doc = pdDocOpen(pdf_path)
        npages = pdDocGetPageCount(doc)
        pdDocClose(doc)
        return npages
    catch e
        @warn "Could not determine page count: $e"
        return 1
    end
end
