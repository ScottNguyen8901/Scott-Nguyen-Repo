from flask import Flask, request, jsonify
import sqlite3
import uuid
import os

app = Flask(__name__)
DB_PATH = 'iot_data.db'

# Create the table if it doesn't exist
def init_db():
    with sqlite3.connect(DB_PATH) as conn:
        conn.execute('''
            CREATE TABLE IF NOT EXISTS records (
                id TEXT PRIMARY KEY,
                json TEXT
            )
        ''')
init_db()

@app.route('/', methods=['POST', 'PUT'])
def create():
    if not request.is_json:
        return jsonify({"error": "Content-Type must be application/json"}), 400

    try:
        data = request.get_json(force=True)
    except Exception as e:
        return jsonify({"error": f"Malformed JSON: {str(e)}"}), 400

    if not isinstance(data, dict):
        return jsonify({"error": "JSON body must be an object"}), 400

    item_id = str(uuid.uuid4())
    with sqlite3.connect(DB_PATH) as conn:
        conn.execute("INSERT INTO records (id, json) VALUES (?, ?)", (item_id, str(data)))
    return jsonify({"id": item_id, "data": data}), 201

@app.route('/<item_id>', methods=['GET'])
def read(item_id):
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.execute("SELECT json FROM records WHERE id = ?", (item_id,))
        row = cursor.fetchone()
        if row:
            return jsonify({"id": item_id, "data": eval(row[0])})
    return jsonify({"error": "Item not found"}), 404

@app.route('/<item_id>', methods=['DELETE'])
def delete(item_id):
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.execute("DELETE FROM records WHERE id = ?", (item_id,))
        if cursor.rowcount:
            return jsonify({"status": "deleted"})
    return jsonify({"error": "Item not found"}), 404

@app.errorhandler(404)
def not_found(e):
    return jsonify({"error": "Invalid route"}), 404

@app.errorhandler(405)
def method_not_allowed(e):
    return jsonify({"error": "Method not allowed"}), 405

@app.errorhandler(500)
def internal_error(e):
    return jsonify({"error": "Internal server error"}), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=80)
