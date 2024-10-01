# AllUNGo: A novel strategy for improving multi-class and multi-label protein function annotation by biological inference and deep learning methods

AllUNGo is a novel strategy designed to improve multi-class and multi-label protein function annotation by integrating biological inference with deep learning methods. Utilizing a model-theoretic approach, AllUNGo constructs ontology embeddings and integrates them with neural networks to predict protein functions. A key feature of AllUNGo is its ability to make zero-shot predictions, allowing the prediction of protein functions for which no proteins in the training set have been annotated. By leveraging formal axioms from the Gene Ontology (GO), AllUNGo significantly enhances annotation accuracy, particularly for GO classes with minimal or no experimental annotations.

This repository contains script which were used to build and train the
AllUNGo model together with the scripts for evaluating the model's
performance.

## Dependencies
* The code was developed and tested using python 3.9.
* To install python dependencies run:
  `pip install -r requirements.txt`
* Install [diamond](https://github.com/bbuchfink/diamond) program on your system (diamond command should be available)
