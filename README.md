# trscc

Run the bootstrap script to perform the requested TRCSS setup flow:

```bash
chmod +x /home/runner/work/trscc/trscc/setup_trcss_repo.sh
/home/runner/work/trscc/trscc/setup_trcss_repo.sh
```

The script will:
1. Check `git` and `gh` availability.
2. Ask you to log in with `gh auth login` if needed.
3. Run `git init` only when the folder is not yet a Git repository.
4. Stage and commit with message `Initial commit: TRCSS v1.0.0-submission`.
5. Ask for confirmation **before any public push**.
6. Run `gh repo create trcss --public --source=. --push`.
7. Tag the commit as `v1.0.0-submission` and push the tag.