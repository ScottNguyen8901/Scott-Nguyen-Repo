import argparse
from rag_chatbot import create_rag, answer_question


QUESTIONS = [
    "What is this homework about?",
    "What are the main steps in building the RAG system?",
    "What is retrieval augmented generation?",
    "What embedding model is suggested in the homework?",
    "What vector database is used?",
    "What are two ways to improve the RAG system?",
    "What percentage is assigned to performance improvement?",
    "What is the capital of France?",
]


def main():
    parser = argparse.ArgumentParser(description="Automated RAG question runner")

    parser.add_argument("--pdf", default="data/ece-595-hw-3.pdf")
    parser.add_argument("--chunk-size", type=int, default=1000)
    parser.add_argument("--chunk-overlap", type=int, default=150)
    parser.add_argument("--k", type=int, default=4)
    parser.add_argument("--llm", default="mistral")
    parser.add_argument("--embedding-model", default="all-MiniLM-L6-v2")
    parser.add_argument("--output", default="rag_results.txt")

    args = parser.parse_args()

    llm, retriever = create_rag(args)
    chat_history = []

    with open(args.output, "w", encoding="utf-8") as f:
        header = f"""
RAG Automated Test Results
==========================

Chunk size: {args.chunk_size}
Chunk overlap: {args.chunk_overlap}
Retriever k: {args.k}
LLM: {args.llm}
Embedding model: {args.embedding_model}

"""
        print(header)
        f.write(header)

        for i, question in enumerate(QUESTIONS, start=1):
            print(f"Running question {i}/{len(QUESTIONS)}")
            print(f"Question {i}: {question}")

            answer, docs = answer_question(llm, retriever, question, chat_history)
            chat_history.append((question, answer))

            console_output = f"""
Answer:
{answer.strip()}

Retrieved source pages:
"""
            print(console_output)
            f.write(f"Question {i}: {question}\n")
            f.write(console_output)

            for j, doc in enumerate(docs, start=1):
                page = doc.metadata.get("page", "unknown")
                source_line = f"  Source {j}: page {page}"
                print(source_line)
                f.write(source_line + "\n")

            separator = "\n" + "-" * 60 + "\n"
            print(separator)
            f.write(separator + "\n")

    print(f"Done. Results saved to {args.output}")


if __name__ == "__main__":
    main()
