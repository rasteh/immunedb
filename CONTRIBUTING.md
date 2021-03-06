# How to contribute
We welcome contributions to ImmuneDB, whether small or large.  For non-developers,
reporting bugs and adding to documentation are both valuable ways to assist.
Developers can help by submitting pull requests for open issues or implementing
new features.

## Reporting Bugs
When reporting a bug make sure to include details about the environment being
used and any error messages.

## Submitting Pull Requests
Pull requests can be submitted for ImmuneDB by anyone.  Before submitting a pull
request, you should make sure they pass local tests by running:

        $ pip install coverage nose
        $ ./tests/run.sh

The environment variable `DB_ADMIN_PASS` should be set to your MySQL
root password and `LL_PATH` should be set to the path of `needleman_wunsch`.

For new features, please add tests as necessary.  To change the regression test
reference database, when running the tests, set the `GENERATE` environment
variable.  Please verify that the changes to the files in `tests/data/regression/`
accurately represent your intentions.

Further, we will **only** accept pull requests if they are
[PEP8](https://www.python.org/dev/peps/pep-0008/) compliant.  To check your code
for compliance, you should get no style warning after running:

        $ pip install pep8
        $ pep8 immunedb

After completing these steps, please submit a pull request.  Tests will be run
on [Travis](https://travis-ci.com/arosenfeld/immunedb).  After they pass, we'll
consider merging your changes or discuss further modifications.
