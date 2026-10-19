# Uploading to GitHub

The empty repository is https://github.com/alvarohumberto/PlasmaWaveLab.

## Browser (drag and drop)

1. Open the repository and click **uploading an existing file**.
2. Drag the **contents** of the `PlasmaWaveLab/` folder (not the folder itself,
   not the ZIP) into the page.
3. **Hidden items are often skipped** by drag and drop: check that
   `.github/workflows/tests.yml`, `.gitignore` and `.gitattributes` appear in
   the list. If not, create `.github/workflows/tests.yml` with
   **Add file → Create new file** and paste its content.
4. Commit to `main`, then open the **Actions** tab and confirm the tests pass.

## Git (alternative)

```bash
git remote add origin https://github.com/alvarohumberto/PlasmaWaveLab.git
git push -u origin main
```

## Release

After CI is green: **Releases → Draft a new release**, tag `v0.3.0`.
For a DOI, enable the repository in Zenodo (GitHub integration) *before*
publishing the release, then add the DOI to `CITATION.cff` and the manuscript.
