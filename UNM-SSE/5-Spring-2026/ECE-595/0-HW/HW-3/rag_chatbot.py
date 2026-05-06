import os
import argparse
from tqdm import tqdm

from langchain_text_splitters import RecursiveCharacterTextSplitter
from langchain_community.document_loaders import PyPDFLoader
from langchain_community.vectorstores import Chroma
from langchain_huggingface import HuggingFaceEmbeddings
from langchain_community.llms import Ollama


PDF_PATH = "data/ece-595-hw-3.pdf"


def build_vector_db(pdf_path, chunk_size, chunk_overlap, embedding_model):
    persist_dir = f"chroma_db_chunk{chunk_size}_overlap{chunk_overlap}"

    print("Loading embedding model...")
    embeddings = HuggingFaceEmbeddings(model_name=embedding_model)

    if os.path.exists(persist_dir):
        print(f"Loading existing vector DB: {persist_dir}")
        return Chroma(
            persist_directory=persist_dir,
            embedding_function=embeddings,
        )

    if not os.path.exists(pdf_path):
        raise FileNotFoundError(f"Could not find PDF at: {pdf_path}")

    print("Loading PDF...")
    pages = PyPDFLoader(pdf_path).load()
    print(f"Loaded {len(pages)} PDF pages.")

    splitter = RecursiveCharacterTextSplitter(
        chunk_size=chunk_size,
        chunk_overlap=chunk_overlap,
    )

    print("Splitting pages into chunks...")
    documents = []
    for page in tqdm(pages, desc="Chunking PDF"):
        documents.extend(splitter.split_documents([page]))

    print(f"Total chunks: {len(documents)}")

    print(f"Creating new vector DB: {persist_dir}")
    vectordb = Chroma.from_documents(
        documents=documents,
        embedding=embeddings,
        persist_directory=persist_dir,
    )

    return vectordb


def answer_question(llm, retriever, question, chat_history):
    docs = retriever.invoke(question)

    context = "\n\n".join(
        f"[Source {i + 1} | Page {doc.metadata.get('page', 'unknown')}]\n{doc.page_content}"
        for i, doc in enumerate(docs)
    )

    history_text = "\n".join(
        f"User: {q}\nAssistant: {a}" for q, a in chat_history[-3:]
    )

    prompt = f"""
You are a precise RAG chatbot.

Use ONLY the document context to answer the question.
If the answer is not clearly in the document context, say:
"I do not know based on the document."

Keep the answer clear and concise.

Conversation history:
{history_text}

Document context:
{context}

Question:
{question}

Answer:
"""

    answer = llm.invoke(prompt)
    return answer, docs


def create_rag(args):
    vectordb = build_vector_db(
        pdf_path=args.pdf,
        chunk_size=args.chunk_size,
        chunk_overlap=args.chunk_overlap,
        embedding_model=args.embedding_model,
    )

    retriever = vectordb.as_retriever(search_kwargs={"k": args.k})
    llm = Ollama(model=args.llm, num_predict=80, temperature=0.2)

    return llm, retriever


def main():
    parser = argparse.ArgumentParser(description="RAG chatbot for ECE 595 HW3")

    parser.add_argument("--pdf", default=PDF_PATH)
    parser.add_argument("--chunk-size", type=int, default=1000)
    parser.add_argument("--chunk-overlap", type=int, default=150)
    parser.add_argument("--k", type=int, default=4)
    parser.add_argument("--llm", default="llama3.2:1b")
    parser.add_argument("--embedding-model", default="all-MiniLM-L6-v2")

    args = parser.parse_args()

    print("RAG Configuration")
    print("-----------------")
    print(f"PDF: {args.pdf}")
    print(f"Chunk size: {args.chunk_size}")
    print(f"Chunk overlap: {args.chunk_overlap}")
    print(f"Retriever k: {args.k}")
    print(f"LLM: {args.llm}")
    print(f"Embedding model: {args.embedding_model}")
    print()

    llm, retriever = create_rag(args)
    chat_history = []

    print("\nRAG chatbot is ready.")
    print("Ask questions about the PDF.")
    print("Type 'exit' to quit.\n")

    while True:
        question = input("Question: ").strip()

        if question.lower() in ["exit", "quit"]:
            print("Goodbye.")
            break

        if not question:
            print("Please enter a question.\n")
            continue

        answer, docs = answer_question(llm, retriever, question, chat_history)
        chat_history.append((question, answer))

        print("\nAnswer:")
        print(answer)

        print("\nRetrieved source pages:")
        for i, doc in enumerate(docs, start=1):
            print(f"  Source {i}: page {doc.metadata.get('page', 'unknown')}")
        print()


if __name__ == "__main__":
    main()
