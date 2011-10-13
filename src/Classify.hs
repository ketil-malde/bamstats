
-- | Classify reads in a BAM-file as innie outie rightie leftie,
--   with stats for distances

module Classify where

import Bio.SamTools.Bam
import Data.List (foldl', intercalate)
import Data.Maybe (isNothing)
import Text.Printf (printf)

data ClassStats = CS { count :: !Int, xsum, x2sum, x3sum, x4sum :: !Double } deriving Show
data Classify = Class { innies, outies, lefties, righties :: !ClassStats } deriving Show

-- | Pretty-print a 'Classify' value.
display :: Classify -> String
display c = unlines $ map (intercalate "\t")
  [ ["Alignment","  count"," mean","stdev"," skew"," kurt"]
  ,  "innies   " : disp1 (innies c)
  ,  "outies   " : disp1 (outies c)
  ,  "lefties  " : disp1 (lefties c)
  ,  "righties " : disp1 (righties c)]
  where disp1 cs = printf "%7d" (count cs) : let ct = fromIntegral (count cs)
                                     in map (printf "%5.1f") [xsum cs/ct, 0, 0, 0]


-- | Default value
cdef :: Classify
cdef = Class c c c c
  where c = CS 0 0 0 0 0

-- | Extract info from alignments
classify :: [Bam1] -> Classify
classify = foldl' class1 cdef 

-- | Update data structure with a single alignment
class1 :: Classify -> Bam1 -> Classify
class1 c0 b 
  | isUnmapped b = c0
  | isOpposite b =
      (if firstUpstream b then add_innie else add_outie) b c0
  | otherwise =
      (if firstUpstream b then add_rightie else add_leftie) b c0

add_innie, add_outie, add_rightie, add_leftie :: Bam1 -> Classify -> Classify
add_innie b c = c { innies = insert b (innies c) }
add_outie b c = c { outies = insert b (outies c) }
add_rightie b c = c { righties = insert b (righties c) }
add_leftie b c = c { lefties = insert b (lefties c) }

insert :: Bam1 -> ClassStats -> ClassStats
insert b (CS c s s2 s3 s4) = CS (c+1) (s+d) (s2+d^2) (s3+d^3) (s4+d^4) 
  where d = maybe (error ("no insert size?\n"++show b)) fromIntegral $ insertSize b

isUnmapped, isOpposite, firstUpstream, isBefore :: Bam1 -> Bool
isUnmapped = isNothing . insertSize -- isUnmap b || isMateUnmap b || mateTargetID b /= targetID b -- <- apparently not the same thing?
isOpposite b = isReverse b `xor` isMateReverse b 
firstUpstream b = isReverse b `xor` isBefore b
isBefore b = position b < matePosition b

xor :: Bool -> Bool -> Bool
x `xor` y = (x || y) && not (x && y) --sigh: xor = (/=)