.PHONY: clean

clean: 
	@python ./scripts/clean.py
	
train: 
	@python ./scripts/train.py
