import json
from http.server import BaseHTTPRequestHandler, HTTPServer

class Handler(BaseHTTPRequestHandler):
    def do_POST(self):
        content_length = int(self.headers.get('Content-Length', 0))
        body = self.rfile.read(content_length).decode()

        if self.path == "/status":
            print(f"\nReceived POST to /status")
            print("----- BEGIN PAYLOAD -----")
            print(body)
            print("------ END PAYLOAD ------")
            self.send_response(200)
            self.end_headers()

        elif self.path == "/program":
            try:
                data = json.loads(body)
                with open("schedule.txt", "w") as f:
                    for entry in data:
                        f.write(f"{entry['time']} {entry['temp']:.1f}\n")

                print(f"\nUpdated schedule.txt from POST to /program")
                print("----- BEGIN SCHEDULE -----")
                with open("schedule.txt") as f:
                    print(f.read().strip())
                print("------ END SCHEDULE ------")

                self.send_response(200)
                self.end_headers()
            except Exception as e:
                self.send_response(400)
                self.end_headers()
                print(f"Error parsing /program POST: {e}")
        else:
            self.send_response(404)
            self.end_headers()

# Start server on localhost:8000
print("HTTP test server running on http://localhost:8000 ...")
HTTPServer(('localhost', 8000), Handler).serve_forever()
