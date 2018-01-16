#!/bin/bash

# ssh-askpass-git.sh uses these as environmental variabls. Its the
# only way I can think of to make that script modular (ie for other accounts)
# since I can't pass arguments to it.
export KW_FOLDER="ssh-key-passwords"
export KW_IDENT="git@rabraker.com"

SSH_KEYFILE="/home/arnold/.ssh/id_rsa"


USAGE="USAGE: sftp-watchdog.sh file-to-watch"

function cleanup {
    if [[ -v sshadd_status ]]; then
        ssh-add -d "$SSH_KEYFILE"
        ssh-agent -k
    fi
}
trap cleanup ERR EXIT

# Make sure an argument is passed
if [[ $# = 0 ]]
then
    echo "$USAGE"
    exit 1
else
    # remove trailing slash
    FILE="$1"
    # echo "$FILE"
fi



upload () {
    FILE="$1"
    # python keyring requires dbus and evidently checks for the environmental
    # variable DISPLAY. We don't need that, set it to zero. This has purportadley
    # been fixed: https://github.com/jaraco/keyring/pull/182
    # but I'm getting the error message anyway.
    export DISPLAY=:0

    SFTP_COMMAND="put ${FILE} /var/www/html/rabraker.com/public_html/matlab-logs/${FILE}"

    # SSH_ASKPASS will only be called if it's not in a terminal. So the only purpose
    # of the pipe is to put it in a bash shell which has not terminal.
    echo `ps -e|grep ssh-agent` > /dev/null
    #echo "before"

    # We have to start our own ssh-agent in cron as far as I know because cron is started
    # OUTSIDE of all sessions and so information about SSH_AUTH_SOCK, e.g., is unknown to
    # to cron.
    eval `ssh-agent -s`

    echo null|SSH_ASKPASS=/home/arnold/local/bin/ssh-askpass-git.sh DISPLAY= ssh-add "$SSH_KEYFILE"

    sshadd_status="$?"

    # We have to run sftp in batch mode so we can get an exit status. So create a temp file,
    # echo our command into it, and delete.
    BATCH_FILE="/tmp/$(basename $)).$$.tmp"

    echo "$SFTP_COMMAND">"$BATCH_FILE"
    STATUS_MESSAGE=`sftp -b "$BATCH_FILE" git@rabraker.com`
    STATUS=$?
    rm "$BATCH_FILE"
    # If sftp fails, send an email and exit so we dont keep failing on the server.
    if [[ "$STATUS" -eq 1 ]]; then
      printf '%s\n' "From: arnold@rabraker.com
To: abraker@fastmail.com
Subject: Watchdog on ${FILE} aborted
Matlab Watchdod aborted because sftp exited with status 1:
${STATUS_MESSAGE}"|sendmail -t -a arnold-rabraker

      exit 1
    fi

    cleanup
}

# upload "$FILE"

LTIME=`stat -c %Z "$FILE"`
while true
do
    ATIME=`stat -c %Z "$FILE"`

    if [[ "$LTIME" != "$ATIME" ]]; then
        echo "uploading..."
        upload "$FILE"
        LTIME="$ATIME"
    fi
    sleep 600
done
