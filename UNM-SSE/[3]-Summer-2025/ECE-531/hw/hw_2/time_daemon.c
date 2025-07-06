// cc -o time_daemon time_daemon.c
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <syslog.h>
#include <time.h>
#include <signal.h>

static volatile int running = 1;

void signal_handler(int sig) {
    if (sig == SIGTERM || sig == SIGINT) {
        running = 0;
    }
}

int main() {
    pid_t pid, sid;

    // Fork the parent process
    pid = fork();
    if (pid < 0) {
        exit(EXIT_FAILURE);
    }
    if (pid > 0) {
        // Exit the parent process
        exit(EXIT_SUCCESS);
    }

    // Create a new session ID for the child process
    sid = setsid();
    if (sid < 0) {
        exit(EXIT_FAILURE);
    }

    // Change working directory to root
    if ((chdir("/")) < 0) {
        exit(EXIT_FAILURE);
    }

    // Close standard file descriptors
    close(STDIN_FILENO);
    close(STDOUT_FILENO);
    close(STDERR_FILENO);

    // Setup signal handling to allow clean shutdown
    signal(SIGTERM, signal_handler);
    signal(SIGINT, signal_handler);

    // Open connection to syslog
    openlog("time_daemon", LOG_PID, LOG_DAEMON);

    while (running) {
        time_t now = time(NULL);
        if (now == ((time_t) -1)) {
            syslog(LOG_ERR, "Failed to get current time");
        } else {
            struct tm *ptm = localtime(&now);
            if (ptm == NULL) {
                syslog(LOG_ERR, "Failed to convert time to local time");
            } else {
                char time_str[64];
                if (strftime(time_str, sizeof(time_str), "%Y-%m-%d %H:%M:%S", ptm) == 0) {
                    syslog(LOG_ERR, "strftime returned 0");
                } else {
                    syslog(LOG_INFO, "Current time: %s", time_str);
                }
            }
        }
        sleep(1);
    }

    syslog(LOG_INFO, "Daemon terminated");
    closelog();

    return EXIT_SUCCESS;
}
