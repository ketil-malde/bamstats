-- | BAM - extract information from BAM files

{-# Language DeriveDataTypeable #-}

module Main where
import Bio.SamTools.Classify
import Bio.SamTools.Bam
import System.Console.CmdArgs

data Options = Stats { numrd :: Maybe Int, inputs :: [FilePath] }
             | Hist { numrd :: Maybe Int, bins, maxdist :: Int, inputs :: [FilePath] }
             | Dump { numrd :: Maybe Int, inputs :: [FilePath] }
             deriving (Data,Typeable,Read,Show)

clfy, hst, dump :: Options
clfy = Stats { numrd = Nothing &= help "max number of reads (default: all)"
                , inputs = []     &= args &= typ "BAM file(s)"
                } &= help "Calculate statistics on insert sizes"
hst  = Hist { numrd = Nothing &= help "max number of reads to include"
            , bins  = 20      &= help "number of bins" &= typ "INT"
            , maxdist = 1000  &= help "max insert size to count" &= typ "INT"
            , inputs = []     &= args &= typ "BAM file(s)"
            } &= help "Collect a histogram of insert sizes"
dump = Dump { numrd = Just 100 &= help "max number of reads (default: 100)" 
            , inputs = []     &= args &= typ "BAM file(s)"            
            } &= help "Dump alignments in the different classes."
       
main :: IO ()
main = do
  o <- cmdArgs $ modes [clfy,hst,dump] 
       &= help "Extract information from BAM files" 
       &= program "bam" &= summary "bam v0.0, ©2011 Ketil Malde"
  let geninp f = (case numrd o of Just x -> take x; Nothing -> id) `fmap` readBams f
      genout = case o of Stats {} -> putStrLn . display . (classify :: [Bam1] -> Stats ClassStats)
                         Hist {}     -> const (print o)
                         Dump {} -> putStrLn . display . (classify :: [Bam1] -> Stats Collect)
  mapM_ (\f -> genout =<< geninp f) $ inputs o