# Contributing to MPlusA
Thank you for taking interest in contributing to the project. This text contains guidance in terms of bug reporting and participating in the development of the library.

## Questions, discussion and enhancements proposals
Any questions related to the library are meant to be issued via Github's issues functionality with the *question* label. Discussion of any topics will be done in subsequent answers to a given issue. The same idea goes for any enhancements proposals: create an issue with the *enhancement* label. The author of the library also accepts direct e-mails; the addess can be found on his personal website or under his profile.

## Bug reporting
When it comes to bug reporting, when creating a Github issue labeled as *bug*, it is advised to include at least (if applicable) the error message and the snippets of code which caused the issue. It is appreciated if the description contains any additional clues regarding the reproduction steps of the problem.

## Contributing code
There are a few suggestions when it comes to contributing code to the library. All programming contributions are appreciated.

### Architecture
The library's architecture is simple and does not need unnecessary complication. The code in `src/mplusa/` directory should be applicable to both tropical and artic algebra, while the code applicable to only one of these should go into `src/mplusa/minplus` or `src/mplusa/maxplus` respectively. This structure may change in the future due to applicability of metaprogramming in tropical-arctic division but currently it should be followed.

When contributing code that is defined only for one of the algebras, it is appreciated if the reason for it is provided.

### Code style
The code style should follow the `PEP 8` standard. This is not a requirement but a suggestion. It is appreciated if the code has docstrings.

### Creating a pull request
The pull request should contain a list with a brief description of changes to the library. The person assigned to it should be the main author of the library. If it's possible, one should provide sources which will allow for easier confirmation of the code's correctness. The branch to which the code will be merged should be the branch named after the future version of the library (currently it will be `0_0_5`) and the source branch should be deleted after the merge.

## Attribution
All people who will contribute a significant new part to the library and will express the desire to be included, will be noted as co-authors in the project's details.

*This guide will get updated with the library's growth.*
