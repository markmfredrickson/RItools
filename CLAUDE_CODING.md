This file contains guidance for Claude when writing  code of any kind, whether
it is for R packages (C++ and R usually), or papers (Python, Bash, R, LaTeX,
usually).

## Setup

Before working on a task, read all files that are relevant
to the task.

For example, in an R package read all files in `R/` and `src/`, `tests/`,
`inst/`, etc.. plus NAMESPACE and DESCRIPTION. When uncertain about scope, read
broadly rather than narrowly --- it is better to understand the surrounding
code than to make changes in isolation. Also read any vignettes that touch the
same functionality. The package uses roxygen2 documentation.

For code written to support the writing of a paper, you will need to understand
*the point* of the code. So that means you need to read the paper files, the
handoff documents, memos, and any other files that will teach you about the
point of the action. When in doubt, ask me about whether your
understanding of the point of the code is correct.

All code written should be able to be saved in files and future researchers on
diverse computing environments should be able to replicate. I'm willing to
require people to use unix/linux at this stage in the game. Some code need not
be permanently tracked in the github repository for a project, other code is
essential for replication.

## Code style

First, prefer boring code over clever code for readability and maintainability.
This does not mean prefer for() loops over other approaches --- since
vectorization is both much faster and a core piece of R language use. In fact,
often avoiding a for() loop can make code more clear than a for() loop because
of all of the overhead necessary to setup the objects to modify, etc.

Second, comment the code explaining **why** the particular sections and lines of
code are there more than **what** they are doing.

## File organization

Group functions by conceptual purpose, one file per coherent unit. Do not let
files grow past roughly 300 lines without checking in with me. When a new
function is conceptually distinct from the existing contents of a file
(different mechanism, different external dependency, different layer of
abstraction), create a new file rather than appending.

## Testing

Write well commented unit tests **before** writing code and/or refactoring. For
example, when working on an R package, write unit tests for tests/testthat
**before** writing code and/or refactoring. Tests should represent both that
the code runs but most importantly they should represent the statistical
principles underlying and justifying the code (if the code is about statistics)
and the **substantive point** of the code being writen. If I am writing code to
square numbers then I should have tests of squaring numbers --- these are more
important tests than tests that when I provide a numeric input I get a numeric
output. Tests should pass before we can judge the code to be correct.

It is not appropriate to remove failing tests rather than fixing source code or
asking for clarification from me. You can skip tests when you can't quickly
resolve failures, but you will need to remember these as future tasks.

Always prioritize readable, maintainable tests over comprehensive coverage.

## Build discipline

When working on an R package, run `devtools::document()` after adding or
modifying roxygen2 documentation. Run `devtools::check()` before considering a
task complete. When adding new exported functions, bump the patch version in
DESCRIPTION (e.g., 0.0.3.0 → 0.0.3.1).

Build discipline will vary for other kinds of code tasks. Many of them will use
a Makefile. So, you should be able to run the makefile for all dependencies for
a given code change.

## Workflow with me

Pause for my review at these checkpoints:
1. After writing tests and before writing implementation.
2. After writing implementation and before running `devtools::check()`.
3. Whenever a design decision arises that the plan does not already resolve.

Do not proceed past a checkpoint without my input.

