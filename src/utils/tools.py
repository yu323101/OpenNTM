"""Tool management utilities for MTB-Evo."""

import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass

from src.utils.logging_config import get_logger

logger = get_logger()


@dataclass
class ToolInfo:
    """Information about a tool."""
    name: str
    path: Optional[Path]
    available: bool
    version: Optional[str] = None


class ToolManager:
    """Manage external bioinformatics tools."""
    
    REQUIRED_TOOLS = ["sickle", "bwa", "samtools", "java"]
    OPTIONAL_TOOLS = ["fastp"]
    
    def __init__(self):
        self._tools: Dict[str, ToolInfo] = {}
        self._check_all()
    
    def _check_all(self):
        """Check all required tools."""
        for tool in self.REQUIRED_TOOLS + self.OPTIONAL_TOOLS:
            path = shutil.which(tool)
            version = None
            
            if path:
                version = self._get_version(tool, path)
                logger.debug(f"Found {tool}: {path} (version: {version or 'unknown'})")
            else:
                logger.warning(f"Tool not found: {tool}")
            
            self._tools[tool] = ToolInfo(
                name=tool,
                path=Path(path) if path else None,
                available=path is not None,
                version=version
            )
    
    def _get_version(self, tool: str, path: str) -> Optional[str]:
        """Get tool version."""
        version_flags = {
            "bwa": "--version",
            "samtools": "--version",
            "java": "-version",
        }
        
        try:
            flag = version_flags.get(tool, "--version")
            result = subprocess.run(
                [path, flag],
                capture_output=True,
                text=True,
                timeout=5
            )
            # 简单解析第一行
            output = result.stdout if result.stdout else result.stderr
            if output:
                first_line = output.strip().split('\n')[0]
                return first_line[:100]  # 限制长度
        except Exception:
            pass
        return None
    
    def check_required(self) -> List[str]:
        """Return list of missing required tools."""
        return [
            name for name, info in self._tools.items()
            if name in self.REQUIRED_TOOLS and not info.available
        ]
    
    def get_path(self, tool: str) -> Path:
        """Get tool path."""
        if tool not in self._tools or not self._tools[tool].available:
            raise RuntimeError(f"Tool not found: {tool}")
        return self._tools[tool].path
    
    def is_available(self, tool: str) -> bool:
        """Check if a tool is available."""
        return tool in self._tools and self._tools[tool].available
    
    def print_status(self):
        """Print tool status to logger."""
        logger.info("Tool Status:")
        for tool in self.REQUIRED_TOOLS + self.OPTIONAL_TOOLS:
            info = self._tools[tool]
            status = "✓" if info.available else "✗"
            path_str = str(info.path) if info.path else "NOT FOUND"
            logger.info(f"  {status} {tool}: {path_str}")
    
    def validate_all(self) -> bool:
        """Validate all required tools are available.
        
        Returns:
            True if all tools available, False otherwise
        """
        missing = self.check_required()
        if missing:
            logger.error(f"Missing required tools: {', '.join(missing)}")
            logger.error("Please activate conda environment: conda activate mtb-evo")
            return False
        
        logger.info("All required tools found")
        return True
