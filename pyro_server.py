import Pyro4
import subprocess
import threading
import os

@Pyro4.expose
class JobServer(object):
    def __init__(self):
        self.lock = threading.Lock()

    def run_job(self, file_name, pdb2pqr):
        def job_thread():
            log_file = "/srv/www/dnaprodb.usc.edu/DNAProDB_v3_frontend/scripts/job_log.txt"
            with open(log_file, "a") as log:
                command = ["/srv/www/dnaprodb.usc.edu/DNAProDB_v3_frontend/scripts/newRunJob.sh", file_name, str(pdb2pqr)]
                subprocess.Popen(command, stdout=log, stderr=log, preexec_fn=os.setpgrp)

        with self.lock:
            threading.Thread(target=job_thread).start()

        return "Job started"

def main():
#    Pyro4.Daemon.serveSimple(
#        {
#            JobServer: "example.jobserver"
#        },
#        ns=False
#    )
    job_server = JobServer()
    daemon = Pyro4.Daemon(port=9090)  # Use port 9090
    uri = daemon.register(job_server, "example.jobserver")
    print("Object <class '__main__.JobServer'>:")
    print("    uri =", uri)
    print("Pyro daemon running.")
    daemon.requestLoop()


if __name__ == "__main__":
    main()
