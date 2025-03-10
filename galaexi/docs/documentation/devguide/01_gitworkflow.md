# Git Workflow

FLEXI uses Git for version management. This section provides an overview of the expectations of a FLEXI developer in regards to the use of Git.
This guide does not endeavor to provide an introduction to Git itself. If you do not know how to use Git, a simple internet search will yield more information
and tutorials than you will ever need.

## IAG GitLab Repository

The source code repository for FLEXI is hosted locally on the IAG GitLab server. The repository can be found [here](https://gitlab.iag.uni-stuttgart.de/flexi/flexi).

## Branches

The `master` branch serves as the main release branch and is protected.

Any new developments must be made in a separate development branch. Development branches are created from the `master` branch with the naming convention `<branch-type>.<branch-name>`, where `<branch-type>` is one of the options in the table below and `<branch-name>` is a one- or two-word description of the changes in the branch.

|  Branch Type   |                                          Description                                                |
| :------------: | :-------------------------------------------------------------------------------------------------: |
|   `feature`    | For large new additions to the code, such as implementing a new physics model.                      |
|  `improvement` | For improvements to existing features, such as adding functionality or refactoring for performance. |
|    `bugfix`    | For fixing a known bug in the code.                                                                 |
|     `test`     | For one-off changes meant for testing. Will typically never be merged into `master`.                |

When choosing a `<branch-type>` from the above table, there may be (many) instances when the classification of a branch falls somewhere between two types. 
In these cases, always choose the higher classification of the two options.

## Merge Requests

Code in development branches can only be merged into `master` via a Merge Request. Merge Requests can be created in the GitLab web interface. Once a Merge Request has been created, it must be assigned at least one other developer for review. The suggested code must receive approval from that reviewer and also pass at least the `Master_Merge` CI/CD Reggie pipeline. To run this pipeline, navigate to Build --> Pipeline Schedules in the FLEXI GitLab repository's web interface. Once there, press the "person" button to take ownership of the pipeline, change the target branch to your development branch and then press the "play" button to run the tests. The status of the testing pipeline can be monitored under Build --> Pipelines.

Note that most large new features will require more scrutiny before merging. Such changes are typically discussed among all developers in a bi-weekly meeting and subjected to more rigorous testing. Often times such developments that implement major new features will also require extension of the testing pipeline itself.
If this is the case, it is **REQUIRED** that you add new unit and reggie tests for your feature. Information adding new tests can be found in the section on [Testing](#09_testing).

## Commits

If commits are to be merged into the `master` branch later, they must be sufficiently documented with a meaningful commit message. All commit messages should be written in English and contain a brief and clear description of the changes made in the commit.

## Issues

TO-DO!

## Milestones and Releases

TO-DO!