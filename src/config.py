"""Configuration management for MTB-Evo."""

import os
from pathlib import Path
from typing import List, Optional

import yaml
from pydantic import BaseModel, Field


class PipelineConfig(BaseModel):
    """Pipeline configuration with command-line priority.
    
    Priority order (highest to lowest):
    1. Command-line arguments
    2. Environment variables
    3. YAML config file
    4. Default values
    """
    
    # 样本和路径
    samples: Path = Field(description="Path to sample list file")
    output_dir: Path = Field(default=Path("results"))
    
    # 线程配置
    threads: int = Field(default=4, ge=1, description="Number of threads for bwa")
    sort_threads: int = Field(default=2, ge=1, description="Number of threads for samtools sort")
    
    # 阈值配置
    min_depth: int = Field(default=10, ge=1, description="Minimum sequencing depth")
    min_coverage: float = Field(default=0.95, ge=0, le=1, description="Minimum genome coverage")
    min_var_freq: float = Field(default=0.75, ge=0, le=1, description="Minimum variant frequency")
    
    # 参考数据路径
    reference_h37rv: Path = Field(default=Path("data/tb_h37rv.fasta"))
    reference_ancestor: Path = Field(default=Path("data/tb.ancestor.fasta"))
    ppe_list: Path = Field(default=Path("data/PPE_INS_loci_Rv.list"))
    
    @classmethod
    def from_args_and_env(
        cls,
        samples: Path,
        output_dir: Optional[Path] = None,
        threads: Optional[int] = None,
        sort_threads: Optional[int] = None,
        config_file: Optional[Path] = None
    ) -> "PipelineConfig":
        """Create config from command-line args and environment.
        
        Args:
            samples: Sample list file (required)
            output_dir: Output directory (optional, default: results)
            threads: Number of threads (optional, default: from env or 4)
            sort_threads: Sort threads (optional, default: from env or 2)
            config_file: YAML config file (optional)
        
        Returns:
            PipelineConfig instance
        """
        # 从YAML加载默认值（如果提供）
        yaml_data = {}
        if config_file and config_file.exists():
            with open(config_file) as f:
                yaml_data = yaml.safe_load(f) or {}
        
        # 构建配置字典，优先级：args > env > yaml > default
        config_dict = {
            "samples": samples,
            "output_dir": output_dir or Path(yaml_data.get("output_dir", "results")),
            "threads": threads or int(os.getenv("MTB_THREADS", yaml_data.get("threads", 4))),
            "sort_threads": sort_threads or int(os.getenv("MTB_SORT_THREADS", yaml_data.get("sort_threads", 2))),
            "min_depth": yaml_data.get("min_depth", 10),
            "min_coverage": yaml_data.get("min_coverage", 0.95),
            "min_var_freq": yaml_data.get("min_var_freq", 0.75),
        }
        
        # 参考数据路径（可从YAML覆盖）
        if "reference" in yaml_data:
            ref = yaml_data["reference"]
            if "h37rv" in ref:
                config_dict["reference_h37rv"] = Path(ref["h37rv"])
            if "ancestor" in ref:
                config_dict["reference_ancestor"] = Path(ref["ancestor"])
        
        if "scripts" in yaml_data and "ppe_list" in yaml_data["scripts"]:
            config_dict["ppe_list"] = Path(yaml_data["scripts"]["ppe_list"])
        
        return cls(**config_dict)
    
    def validate_paths(self) -> List[str]:
        """Validate that required paths exist.
        
        Returns:
            List of error messages (empty if all valid)
        """
        errors = []
        
        if not self.samples.exists():
            errors.append(f"Sample list not found: {self.samples}")
        
        if not self.reference_h37rv.exists():
            errors.append(f"Reference genome not found: {self.reference_h37rv}")
        
        if not self.reference_ancestor.exists():
            errors.append(f"Ancestor sequence not found: {self.reference_ancestor}")
        
        if not self.ppe_list.exists():
            errors.append(f"PPE list not found: {self.ppe_list}")
        
        return errors
