// cc -D_FORTIFY_SOURCE=2 -O2 -Wall -Wextra -o time_daemon time_daemon.c
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <syslog.h>
#include <time.h>
#include <signal.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/resource.h>

static volatile sig_atomic_t running = 1;

void signal_handler(int sig) {
    if (sig == SIGTERM || sig == SIGINT) {
        running = 0;
    }
}

// Minimal privilege drop function
void drop_privileges() {
    if (setuid(65534) != 0) {  // 65534 = nobody user on most systems
        syslog(LOG_ERR, "Failed to drop privileges: %s", strerror(errno));
        exit(EXIT_FAILURE);
    }
}

// Set resource limits (defensive)
void set_limits() {
    struct rlimit rl;

    rl.rlim_cur = 10 * 1024 * 1024;  // 10 MB max memory usage
    rl.rlim_max = 10 * 1024 * 1024;
    setrlimit(RLIMIT_AS, &rl);

    rl.rlim_cur = 1024;  // max open files
    rl.rlim_max = 1024;
    setrlimit(RLIMIT_NOFILE, &rl);
}

int main() {
    pid_t pid, sid;

    // Fork the parent process
    pid = fork();
    if (pid < 0) {
        perror("fork failed");
        exit(EXIT_FAILURE);
    }
    if (pid > 0) {
        exit(EXIT_SUCCESS);
    }

    // Create a new session ID
    sid = setsid();
    if (sid < 0) {
        perror("setsid failed");
        exit(EXIT_FAILURE);
    }

    // Change working directory to /
    if (chdir("/") < 0) {
        perror("chdir failed");
        exit(EXIT_FAILURE);
    }

    // Close file descriptors
    close(STDIN_FILENO);
    close(STDOUT_FILENO);
    close(STDERR_FILENO);

    // Set restrictive umask
    umask(027);

    // Signal handling
    struct sigaction sa;
    memset(&sa, 0, sizeof(sa));
    sa.sa_handler = signal_handler;
    sigaction(SIGTERM, &sa, NULL);
    sigaction(SIGINT, &sa, NULL);

    // Set resource limits
    set_limits();

    // Open syslog
    openlog("time_daemon", LOG_PID | LOG_CONS, LOG_DAEMON);
    syslog(LOG_INFO, "Daemon started securely");

    // Drop privileges
    drop_privileges();

    while (running) {
        time_t now = time(NULL);
        if (now == ((time_t) -1)) {
            syslog(LOG_ERR, "Failed to get current time");
        } else {
            struct tm ptm;
            char time_str[64];

            if (localtime_r(&now, &ptm) == NULL) {
                syslog(LOG_ERR, "Failed to convert time");
            } else if (strftime(time_str, sizeof(time_str), "%Y-%m-%d %H:%M:%S", &ptm) == 0) {
                syslog(LOG_ERR, "strftime returned 0");
            } else {
                syslog(LOG_INFO, "Current time: %s", time_str);
            }
        }
        sleep(1);
    }

    syslog(LOG_INFO, "Daemon terminated");
    closelog();

    return EXIT_SUCCESS;
}