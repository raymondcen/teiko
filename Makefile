.PHONY: setup pipeline dashboard

setup:
	pip install -r requirements.txt

pipeline:
	python load_data.py

dashboard:
	python -m streamlit run dashboard.py --server.address 0.0.0.0 --server.port 8501
