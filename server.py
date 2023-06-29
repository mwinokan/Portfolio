#!/usr/bin/env python3

import http.server
import socketserver

PORT = 8000

class MyHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_my_headers()
        http.server.SimpleHTTPRequestHandler.end_headers(self)

    def send_my_headers(self):
        self.send_header("Cache-Control", "no-cache, no-store, must-revalidate")
        self.send_header("Pragma", "no-cache")
        self.send_header("Expires", "0")

with socketserver.TCPServer(("", PORT), MyHTTPRequestHandler) as httpd:
    print("Server started at localhost:" + str(PORT))
    httpd.serve_forever()

# class NoCacheHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):

#     def send_response_only(self, code, message=None):
#         super().send_response_only(code, message)
#         self.send_header('Cache-Control', 'no-store, must-revalidate')
#         self.send_header('Expires', '0')

if __name__ == '__main__':
    http.server.test(HandlerClass=MyHTTPRequestHandler,port=PORT)
    # http.server.test(HandlerClass=NoCacheHTTPRequestHandler,port=PORT)
