# SeqMgr
A standalone sequence reading library.

This library provides the functionality about reading sequences information from database.
No dependency is required and the reading procedure is fast. Basically it is implemented in
a way of in-place parsing, so the manager is required to be kept in memory until you don't
need it. The intermediate data structure is simple, which means that you're recommended to
build your own data structure in your program. Comparing with the parsing, the data transfter
should be less costly.

## License
BSD License
