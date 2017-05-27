# Please specify the location of zmdp's libraries and include files
ZMDP_LIB_DIR = /Users/mthon/dev/zmdp/lib/darwin16
ZMDP_INCLUDE_DIR = /Users/mthon/dev/zmdp/include

LIBS = -L$(ZMDP_LIB_DIR) -lzmdpPomdpBounds -lzmdpPomdpCore -lzmdpPomdpParser -lzmdpCommon
INCLUDE = -I$(ZMDP_INCLUDE_DIR)

POMDPS = tiger.pomdp paint.pomdp bridge.pomdp network.pomdp shuttle.pomdp maze4x3.pomdp cheese.pomdp

pomdps: pomdp_to_json
	@python get_original_pomdps.py
	@for p in $(POMDPS); do ./pomdp_to_json $$p > $$p.json; done

pomdp_to_json: pomdp_to_json.cpp
	$(CXX) -O2 $(INCLUDE) pomdp_to_json.cpp -o pomdp_to_json $(LIBS)
