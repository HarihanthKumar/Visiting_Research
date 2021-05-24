#Interval class is to store interval.
#The index attribute is used to the rowNumber of the entry in the 
#methylation database that can be used to access the methylation value
class Interval:
    #initializes the data members.
    def __init__(self,start:int,end:int, index:int):
        self.start=start
        self.end=end
        self.index=index

#IntervalTree class is used to build the tree of interval nodes. It is a Binary Search Tree. 
class IntervalTree:
    #initializes the data members.
    #left and right are pointers to the left and right child interval subtrees. 
    #height of the tree gives the depth.
    #maxValue denotes the maximum between the maximum end value of the left and the right subtrees. This is used to make the decision
    #of which subtree to search. 
    def __init__(self):
        self.left=None
        self.right=None
        self.node=None
        self.maxValue=None
        self.height = None

    #creates a tree node. Initially when there is only one node, the maximum value is the end value of the node itself and the height 
    #is 1.
    def createNode(self,interval:Interval):
        returnNode=IntervalTree()
        returnNode.node=interval
        returnNode.maxValue=interval.end
        returnNode.height = 1
        return returnNode

    #Returns the height of that node in the tree. Height is defined as the position of the node from the bottom of the tree.     
    def getHeight(self, root) -> int:
        if root is None: 
            return 0
        else:
            return root.height

    #returns the balance factor. This is used to balance the tree if tree becomes skewed. Balancing ensures that the height of the 
    #entire tree and all its subtrees is always log(n) where n is the number of nodes in the tree. Balance factor of a node is 
    #defined as the difference between the height of its left subtree and right subtree. 
    def getBalance(self, root) -> int:
        if root is None:
            return 0
        else:
            return (self.getHeight(root.left) - self.getHeight(root.right))
    
    #The two rotate functions are here to ensure balance of a node. These functions are appropriately called when the absolute height 
    #difference between the left and right subtrees of a node is greater than 1 (So, the balance factor can be only one of the three 
    #values [-1,0,1]). These functions rotate the tree about a node. When the trees are rotated the maxValue is updated to make sure 
    #it remains intact and correct. 
    def rotateLeft(self, root):
        rightOfRoot = root.right
        leftOfRight = rightOfRoot.left

        rightOfRoot.left = root
        root.right = leftOfRight
        root.maxValue = max(root.node.end, leftOfRight.maxValue if leftOfRight is not None else root.node.end)

        root.height = (1 + max(self.getHeight(root.left), self.getHeight(root.right)))
        rightOfRoot.height = (1 + max(self.getHeight(rightOfRoot.left), self.getHeight(rightOfRoot.right)))

        return rightOfRoot

    def rotateRight(self, root):
        leftOfRoot = root.left
        rightOfLeft = leftOfRoot.right

        leftOfRoot.right = root
        root.left = rightOfLeft

        leftOfRoot.maxValue = max(leftOfRoot.node.end, root.maxValue if root is not None else leftOfRoot.node.end)

        root.height = (1 + max(self.getHeight(root.left), self.getHeight(root.right)))
        leftOfRoot.height = (1 + max(self.getHeight(leftOfRoot.left), self.getHeight(leftOfRoot.right)))

        return leftOfRoot

    #Used to insert an interval to the node tree. Let the interval to be inserted be [a,b]. The search for insertion starts at the 
    #root and let the interval of the root be [x,y]. If a < x, then the function is recursively called to insert the current interval 
    #somewhere in the left subtree. If the a > x, then then the current interval is inserted somewhere in the right subtree. Inserting 
    # in this careful manner helps immensely when an interval has to be searched because you can make the decision in which direction to
    # search because we already know elements that have smaller start value are to the left and greater start value are to the right. 
    # Since the search space is reduced to half in each step, it only takes log(n) steps to find an element instead of 'n' in case of a
    # linear search. This type of search is called binary search, hence the tree goes by that name. The tree is also checked for balance
    #after each insertion and is rotated to make sure of the validity of the balance factor condition. 
    def insertNode(self,root,interval:Interval):
        if root is None:
            return self.createNode(interval)
        elif(root.node.start>interval.start):
            root.left=self.insertNode(root.left,interval)
        elif(root.node.start<interval.start):
            root.right=self.insertNode(root.right,interval)
        if(root.maxValue<interval.end):
            root.maxValue = interval.end
        root.height = (1 + max(self.getHeight(root.left), self.getHeight(root.right)))
        balance = self.getBalance(root)
        if(balance > 1 and interval.start < root.left.node.start):
            return self.rotateRight(root)
        if(balance < -1 and interval.start > root.right.node.start):
            return self.rotateLeft(root)
        if(balance > 1 and interval.start > root.left.node.start):
            root.left = self.rotateLeft(root.left)
            return self.rotateRight(root)
        if(balance < -1 and interval.start < root.right.node.start):
            root.right = self.rotateRight(root.right)
            return self.rotateLeft(root)
        return root

    #This function is used to check if the two intervals overlap. 
    def overlapCheck(self,interval1:Interval,interval2:Interval)->bool:
        if((interval2.start in range(interval1.start, (interval1.end+1))) or (interval2.end in range(interval1.start, (interval1.end+1)))):
            return True
        return False

    #This function is used to search to find the overlapping intervals. The above overlapping function is used to check if two intervals
    #overlap(The methylation intervals and the mutation interval passed from the for loop). If it overlaps, add the rowNumber of the
    #entry in the mutation database to the list corresponding to that methylation interval. 
    def intervalSearch(self,root,interval:Interval,mutationsList: list, rowNumber: int):
        if root is None:
            return
        if(self.overlapCheck(root.node,interval) is True):
            index = root.node.index
            mutationsList[index].append(rowNumber); 
        if((root.left is not None) and (root.left.maxValue >= interval.start)):
            self.intervalSearch(root.left,interval,mutationsList, rowNumber)
        if((root.right is not None) and (root.right.maxValue >= interval.start)):
            self.intervalSearch(root.right,interval,mutationsList, rowNumber)


#Interval Trees using Binary Search Trees are great and if it stills piques your curiosity, you can visit this great article on them
#to know about them in detail https://www.geeksforgeeks.org/interval-tree/.
