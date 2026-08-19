#!/usr/bin/env bash
set -euo pipefail

COMMIT_MESSAGE="Initial commit: TRCSS v1.0.0-submission"
TAG_NAME="v1.0.0-submission"
REPO_NAME="trcss"

check_tool() {
  local tool_name="$1"
  if ! command -v "${tool_name}" >/dev/null 2>&1; then
    return 1
  fi
}

if ! check_tool git; then
  echo "git is not installed."
  echo "Install it first, then re-run this script:"
  echo "  macOS: brew install git"
  echo "  Ubuntu/Debian: sudo apt-get update && sudo apt-get install -y git"
  exit 1
fi

if ! check_tool gh; then
  echo "GitHub CLI (gh) is not installed."
  echo "Install and authenticate, then re-run this script:"
  echo "  macOS: brew install gh"
  echo "  Ubuntu/Debian: (https://github.com/cli/cli/blob/trunk/docs/install_linux.md)"
  echo "  Then login: gh auth login"
  exit 1
fi

if ! gh auth status >/dev/null 2>&1; then
  echo "gh is installed but you are not logged in."
  echo "Run:"
  echo "  gh auth login"
  exit 1
fi

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "Initializing git repository..."
  git init
else
  echo "Git repository already initialized."
fi

echo "Staging files..."
git add .

if git diff --cached --quiet; then
  echo "No staged changes to commit."
else
  echo "Creating commit: ${COMMIT_MESSAGE}"
  git commit -m "${COMMIT_MESSAGE}"
fi

echo
read -r -p "This will create a PUBLIC GitHub repo and push code/tag. Continue? [y/N] " confirm
if [[ ! "${confirm}" =~ ^[Yy]$ ]]; then
  echo "Aborted before public push."
  exit 0
fi

echo "Creating and pushing public repository ${REPO_NAME}..."
gh repo create "${REPO_NAME}" --public --source=. --push

if git rev-parse "${TAG_NAME}" >/dev/null 2>&1; then
  echo "Tag ${TAG_NAME} already exists locally."
else
  echo "Creating tag ${TAG_NAME}..."
  git tag "${TAG_NAME}"
fi

echo "Pushing tag ${TAG_NAME}..."
git push origin "${TAG_NAME}"

echo "Done."
