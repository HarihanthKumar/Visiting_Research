class Interval:
    def __init__(self,start:int,end:int, index:int):
        self.start=start
        self.end=end
        self.index=index

class IntervalTree:
    def __init__(self):
        self.left=None
        self.right=None
        self.node=None
        self.maxValue=None
        self.height = None

    def createNode(self,interval:Interval):
        returnNode=IntervalTree()
        returnNode.node=interval
        returnNode.maxValue=interval.end
        returnNode.height = 1
        return returnNode
        
    def getHeight(self, root) -> int:
        if root is None: 
            return 0
        else:
            return root.height

    def getBalance(self, root) -> int:
        if root is None:
            return 0
        else:
            return (self.getHeight(root.left) - self.getHeight(root.right))
    
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

    def overlapCheck(self,interval1:Interval,interval2:Interval)->bool:
        if((interval2.start in range(interval1.start, (interval1.end+1))) or (interval2.end in range(interval1.start, (interval1.end+1)))):
            return True
        return False

    def intervalSearch(self,root,interval:Interval,methylationValuesIndices: list):
        if root is None:
            return
        if(self.overlapCheck(root.node,interval) is True):
            rowNumber = root.node.index
            methylationValuesIndices.append(rowNumber)
        if((root.left is not None) and (root.left.maxValue >= interval.start)):
            self.intervalSearch(root.left,interval,methylationValuesIndices)
        if((root.right is not None) and (root.right.maxValue >= interval.start)):
            self.intervalSearch(root.right,interval,methylationValuesIndices)


