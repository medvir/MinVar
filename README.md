## Minority Variants

OptimAssembly was merged as subtree following the instructions
[here](https://help.github.com/articles/about-git-subtree-merges).

Rejected subtree pushes were made possible by nesting git commands as
explained [here](https://coderwall.com/p/ssxp5q). Example

    git push ReportDRM `git subtree push --prefix ReportDRM master`
