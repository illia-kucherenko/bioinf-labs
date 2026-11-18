# Assignment: DNA Pattern Generation from Gene Promoters

## Objective
Implement a software application that generates DNA patterns from DNA sequences in gene promoters by computing the $C+G\%$ and the information content values from each sliding window. To ensure the algorithm implementation works correctly, follow these steps:

## Steps
1. **Use the sequence:**
   $$S = \text{CGGACTGATCTATCT(A27)CGTAGCATCTATCGATCTATCTAGCGATCTA}$$
   *(Note: `(A27)` indicates the nucleotide `A` repeated 27 times, and `CTAGCGATCTA` is the continuation of the sequence written on the line below).*

2. **Use a sliding window:**
   $$L = 30\text{ b}$$

3. **Build a function to process the CpG ($C+G$) content:**
   * Expected return value: $CG \approx 29.27$

4. **Build a function to process the information content:**
   * Expected return value: $I_C \approx 27.53$

5. **Plot the pattern.**

6. **Calculate the center of weight of the pattern.**

7. **Take the center of each pattern and plot it on a second chart.**

8. **Take the DNA sequences of promoters and generate patterns:**
   * Open **PromKappa** from GitHub and check if the generated pattern is the same.

---

## Formulas

For a sequence $S = \{x_1, \dots, x_{|S|}\}$:

$$\sigma(S) = 100 - \left( \frac{\sum_{u=1}^{|S|-1} \left( \frac{\sum_{i=1}^{|S|-u} f(x_i, x_{u+i})}{(|S|-u) \times 100} \right)}{|S|-1} \right)$$

Where:

$$f(x_i, x_{u+i}) = \begin{cases} +1, & x_i = x_{u+i} \\ 0, & x_i \neq x_{u+i} \end{cases}$$