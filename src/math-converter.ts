import { mathjax } from "mathjax-full/js/mathjax";
import { TeX } from "mathjax-full/js/input/tex";
import { SVG } from "mathjax-full/js/output/svg";
import { liteAdaptor } from "mathjax-full/js/adaptors/liteAdaptor";
import { RegisterHTMLHandler } from "mathjax-full/js/handlers/html";
import { AllPackages } from "mathjax-full/js/input/tex/AllPackages";

/**
 * Handles conversion of LaTeX math expressions to SVG.
 * Used for rendering mathematical formulas in VASP Wiki documentation within hover tooltips.
 */
export class MathConverter {
    /**
     * Converts a LaTeX string to an SVG string.
     * Uses MathJax's lite DOM adaptor for server-side rendering.
     * 
     * @param latex - The LaTeX math string (e.g., "E = mc^2").
     * @returns The rendered SVG markup as a string.
     */
    convert(latex: string): string {
        const adaptor = liteAdaptor();
        RegisterHTMLHandler(adaptor);

        // Configure MathJax with TeX input and SVG output
        const html = mathjax.document("", {
            InputJax: new TeX({ packages: AllPackages }),
            OutputJax: new SVG({ fontCache: "none" })
        });

        // Convert the math
        const svgNode = html.convert(latex, { display: false });

        // Serialize the SVG node to string
        return adaptor.innerHTML(svgNode);
    }
}
