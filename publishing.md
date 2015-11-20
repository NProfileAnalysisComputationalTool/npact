# Publishing

Most of the publishing is accomplished by `publish.sh`. It wraps a
`git pull` and does a few extra build steps.

## Server Requirements

In addition to all the requirements listed in the README:

* A proxy server to sit in front of the site
* The included `apache.conf` is a pretty good template
* a `.htpasswd` file to protect `/npact/management`


## First time
The first time NPACT is published to a new server a few more steps are
going to be needed.

* Ensure `ReleaseNotes.md` is updated and the version at the top has todays date of publish
* `git tag v${VERSION}`
* `git push --tags`
* git clone the repository to the desired publish destination
    * `git config core.sharedRepository group`
* Include the `apache.conf` under the desired domain config.
* From development host check variables at the top of `publish.sh` and run `publish.sh`

## Republish
* Ensure `ReleaseNotes.md` is updated and the version at the top has todays date of publish
* `git tag v${VERSION}`
* `git push --tags`
* From development host check variables at the top of `publish.sh` and run `publish.sh`
* Go to http://$HOST/npact/management and make sure the tqdaemon gets restarted
