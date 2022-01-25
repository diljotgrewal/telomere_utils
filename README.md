
# telomere utils  
  
### Search Space for telomeres  
  
  
Grid - Search is done with the following parameters:  
  
  
|Parameter|Setting 1|Setting 2|  
| -- | -- | -- |  
| Search Direction| Start -> End| End -> start|
| Kmers | ['TTAGGG', 'TCAGGG', 'TGAGGG', 'TTGGGG']| ['CCCTAA', 'CCCTGA', 'CCCTCA', 'CCCCAA'] (reverse complemented)|
| Read direction | Left -> Right | Right -> Left|


This gives us 2^3 combinations:

|        k-mers        | Read Direction | Search Direction |
|:--------------------:|:--------------:|:----------------:|
| Regular              | Left -> Right  | Start -> End     |
| Regular              | Right -> Left  | Start -> End     |
| Reverse Complemented | Left -> Right  | Start -> End     |
| Reverse Complemented | Right -> Left  | Start -> End     |
| Regular              | Left -> Right  | End -> Start     |
| Regular              | Right -> Left  | End -> Start     |
| Reverse Complemented | Left -> Right  | End -> Start     |
| Reverse Complemented | Right -> Left  | End -> Start     |


These combinations are represented in the output file as a combination of flags:
|reverse_complemented_kmers|reversed_read|match_beginning_of_read|
|:--------------------:|:--------------:|:----------------:|
|F|F|T|
|F|T|T|
|T|F|T|
|T|T|T|
|F|F|F|
|F|T|F|
|T|F|F|
|T|T|F|


Examples:


Our read sequence:
```
CCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCT
```

the read matches reverse complemented k-mer `CCCCAA`  from beginning to end except for 4 bases at the very end.

#### Use case 1:
*No Telomere Found*

#### Use case 2:
*No Telomere Found*

#### Use case 3:
Since search starts from the beginning, start is asssumed to be 0

First Telomere: [0,5]
Last Telomere: [90,95]

start and end positions of telomere: **0,96**

#### Use case 4:
First we reverse the read:
```
TCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCC
```
Since search starts from the beginning, start is always 0
First Telomere: [6,11]
Last Telomere: [90,95]

reported start and end positions of telomere: **0,96**

#### Use case 5:
*No Telomere Found*

#### Use case 6:
*No Telomere Found*

#### Use case 7:
We're going to match from end of read to start. 

End is always 100  in this case (end of read)

The first telomere we see is [90,95]
The last telomere we see is [0,5]

so start will be set to 0 ( start of last telomere seen)
reported start and end positions of telomere: **0,100**


#### Use case 8:
First we reverse the read:
```
TCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCC
```
End is always 100  in this case (end of read)

The first telomere we see is [90,95]
The last telomere we see is [6,11]

reported start and end positions of telomere: **6,100**



Final output:
```
read_id,sample_id,strand,chromosome,start,end,telomere_start,telomere_end,reverse_complemented_kmers,reversed_read,match_beginning_of_read,readend
RID,SA,-,1,10255,10333,0,96,True,False,True,2
RID,SA,-,1,10255,10333,0,96,True,True,True,2
RID,SA,-,1,10255,10333,0,100,True,False,False,2
RID,SA,-,1,10255,10333,6,100,True,True,False,2
```
